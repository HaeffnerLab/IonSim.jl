using SparseArrays: rowvals, nzrange, spzeros, spdiagm, findnz, rowvals, nonzeros
using FunctionWrappers: FunctionWrapper
using PolynomialRoots: roots
using QuantumOptics: SparseOperator, embed

export hamiltonian

"""
    hamiltonian(
            chamber::Chamber; timescale::Real=1, lamb_dicke_order::Union{Vector{Int},Int}=1,
            rwa_cutoff::Real=Inf, displacement="truncated", time_dependent_eta=false
        )
Constructs the Hamiltonian for `chamber` as a function of time. Return type is a function
`h(t::Real, ψ)` that, itself, returns a `QuantumOptics.SparseOperator`.

**args**
* `timescale`: e.g. a value of 1e-6 will take time to be in ``\\mu s``
* `lamb_dicke_order`: Only consider terms that change the phonon number by up to this value.
    If this is an `Int`, then the cutoff is applied to all modes. If this is a `Vector{Int}`,
    then `lamb_dicke_order[i]` is applied to the iᵗʰ mode, according to the order in
    `basis(chamber)`.
    Note: this isn't quite the same thing as the Lamb-Dicke approximation since setting
    `lamb_dicke_order=1` will retain, for example, terms proportional to ``a^\\dagger a ``.
* `rwa_cutoff`: drop terms in the Hamiltonian that oscillate faster than this cutoff.
* `displacement`: This can be either `"truncated"`(default) or `"analytic"`.

   When an atom is irradiated, both the atom's energy and its momentum will generally be
   affected. For an atom in a harmonic potential, the exchange of momentum can be modeled as
   a displacement operation ``D(α=iηe^{-iνt}) = exp[αa^† - α^*a]``, where ``η`` is the
   Lamb-Dicke parameter, which can be described equivalently as either being proportional to
   the square root of the ratio of the recoil frequency with the ground state energy of the
   atom's motion or as the ratio of the spread of the ground state wavefunction to the
   wavelength of the laser.

   When `"truncated"` is selected, the matrix elements of ``D(α)`` are computed by
   constructing ``α^* a, αa^†`` in a truncated basis (according to the dimension specified in
   your model) and then exponentiating their difference. This has the advantage, amongst
   other things, of guaranting unitarity.

   If `"analytic"` is selected, then the matrix elements are computed assuming an infinite-
   dimensional Hilbert space.

   For small displacements (``η ≪ N``, where ``N`` is the dimension of the motion's Hilbert
   space), both of these methods will be good approximations.
* `time_dependent_eta::Bool`: In addition to impacting the vibrational subspace directly, a
   change in the trap frequency, ``δν``, will also change the Lamb-Dicke parameter. Since
   typically ``δν≪ν``, this effect will be small ``η ≈ η₀(1 + δν/2ν)`` and doesn't warrant
   the additional computational resources needed to calculate and update it in time. In this
   case, we can set `time_dependent_eta=false` (default), which will set ``η(t) = η₀``.

"""
function hamiltonian(
    chamber::Chamber;
    timescale::Real=1,
    lamb_dicke_order::Union{Vector{Int}, Int}=1,
    rwa_cutoff::Real=Inf,
    displacement::String="truncated",
    time_dependent_eta::Bool=false
)
    return hamiltonian(
        chamber,
        iontrap(chamber),
        timescale,
        lamb_dicke_order,
        rwa_cutoff,
        displacement,
        time_dependent_eta
    )
end

#############################################################################################
# Hamiltonian for a linear configuration of ions
#############################################################################################

# At each time step, this function updates in-place the 2D array describing the full system
# Hamiltonian.
function hamiltonian(
    chamber::Chamber,
    iontrap::LinearChain,
    timescale::Real,
    lamb_dicke_order::Union{Vector{Int}, Int},
    rwa_cutoff::Real,
    displacement::String,
    time_dependent_eta::Bool
)
    b, indxs, cindxs = _setup_base_hamiltonian(
        chamber,
        timescale,
        lamb_dicke_order,
        rwa_cutoff,
        displacement,
        time_dependent_eta
    )
    aui, gbi, gbs, bfunc, δνi, δνfuncs = _setup_fluctuation_hamiltonian(chamber, timescale)
    S = SparseOperator(basis(chamber))
    function f(t, ψ)  # a two argument function is required in the QuantumOptics solvers
        @inbounds begin
            @simd for i in 1:length(indxs)
                bt_i, conj_bt_i = b[i](t)::Tuple{ComplexF64, ComplexF64}
                @simd for j in 1:length(indxs[i])
                    i1, i2 = indxs[i][j]
                    S.data[i1, i2] = bt_i
                    if i1 != i2
                        S.data[i2, i1] = conj(bt_i)
                        if length(cindxs[i]) != 0
                            flag = cindxs[i][1][1]
                            i3, i4 = cindxs[i][j+1]
                            if flag == -1
                                S.data[i3, i4] = -conj_bt_i
                                S.data[i4, i3] = -conj(conj_bt_i)
                            else
                                S.data[i3, i4] = conj_bt_i
                                S.data[i4, i3] = conj(conj_bt_i)
                            end
                        end
                    end
                end
            end
            if length(gbi) == 0 && length(δνi) == 0
                return S
            else
                @simd for indx in aui
                    S.data[indx, indx] = complex(0.0)
                end
                @simd for i in 1:length(gbi)
                    zeeman_t = bfunc(t)::Float64
                    @simd for j in 1:length(gbi[i])
                        indx = gbi[i][j]
                        S.data[indx, indx] += zeeman_t * gbs[i][j]
                    end
                end
                @simd for i in 1:length(δνi)
                    δν_t = δνfuncs[i](t)::Float64
                    @simd for n in 1:length(δνi[i])
                        @simd for indx in δνi[i][n]
                            S.data[indx, indx] += n * δν_t
                        end
                    end
                end
            end
        end
        return S
    end
    return f
end

#=
The purpose of the hamiltonian function is to evaluate a vector of time-dependent functions
and use the returned values to update, in-place, a pre-allocated array.

The pre-allocated array holds the full Hamiltonian -- a tensor product defined over all of
the individual ion and vibrational mode subspaces -- at a particular point in time.

However, we don't know a priori the exact form of the Hilbert space or the details of the
Hamiltonian's time dependence, since it will be defined by the user.

The _setup_hamiltonian function extracts this user-defined information from a <:Chamber struct
and converts it into two vector of vectors of indices, corresponding to redundant (see note
below) matrix elements of the Hamiltonian, and a matched vector of time-dependent functions
for updating these elements.

------- Note -------
Since the terms of the Hamiltonian will always be of the form of a single ion operator
tensored with a single vibrational mode operator, there will be a lot of redundancy in the
Hamiltonian's matrix elements. E.g.

                                     [ σ₊ ⊗ D(α(t))      0         ]
             H = 𝐼 ⊗ σ₊ ⊗ D(α(t)) =  [       0        σ₊ ⊗ D(α(t)) ]

So to avoid unnecessarily evaluating functions more than once, _setup_hamiltonian also
returns a vector of vectors of indices that keep track of this redundancy.

Also, we have: <m|D(α)|n> = (-1)^(n-m) × conjugate(<n|D(α)|m>). We keep track of this in an
additional vector of vectors of indices.

Finally, since we require the Hamiltonian to be Hermitian, h[i, j] = conj(h[j, i]), this
function does not keeps track of only one of these pairs.
=#
function _setup_base_hamiltonian(
    chamber,
    timescale,
    lamb_dicke_order,
    rwa_cutoff,
    displacement,
    time_dependent_eta
)
    rwa_cutoff *= timescale
    allmodes = reverse(modes(chamber))
    L = length(allmodes)
    νlist = Tuple([frequency(mode) for mode in allmodes])
    mode_dims = [modecutoff(mode) + 1 for mode in allmodes]

    all_ions = reverse(ions(chamber))
    Q = prod([shape(ion)[1] for ion in all_ions])
    ion_arrays = [spdiagm(0 => [true for _ in 1:shape(ion)[1]]) for ion in all_ions]

    ηm, Δm, Ωm =
        _ηmatrix(chamber), _Δmatrix(chamber, timescale), _Ωmatrix(chamber, timescale)
    lamb_dicke_order = _check_lamb_dicke_order(lamb_dicke_order, L)
    ld_array, rows, vals = _ld_array(mode_dims, lamb_dicke_order, νlist, timescale)
    if displacement == "truncated" && time_dependent_eta
        rootlist = map(x -> real.(roots(_He(x))), mode_dims)
    end

    indxs_dict = Dict()
    repeated_indices = Vector{Vector{Tuple{Int64, Int64}}}(undef, 0)
    conj_repeated_indices = Vector{Vector{Tuple{Int64, Int64}}}(undef, 0)
    functions = FunctionWrapper[]
    work_eta = zeros(Float64, L)
    local ts, ion_rows, ion_cols, ion_idxs, ion_reps, rn

    # iterate over ions and lasers
    for n in eachindex(all_ions), m in eachindex(lasers(chamber))
        if m ≡ 1
            rn = length(all_ions) - n + 1
            ts = subleveltransitions(all_ions[n])
            C = sum([
                i * real.(sigma(all_ions[n], reverse(ts[i])...).data) for i in 1:length(ts)
            ])
            if length(all_ions) == 1
                K = C
            else
                K = kron(ion_arrays[1:(n-1)]..., C, ion_arrays[(n+1):length(all_ions)]...)
            end
            ion_rows, ion_cols, ion_vals = findnz(K)
            ion_idxs = sortperm(real.(ion_vals))
            ion_reps = Int(Q / size(C, 1))
        end
        ηnm = view(ηm, rn, m, :)
        function ηlist(t)
            for i in 1:L
                η = ηnm[i]
                typeof(η) <: Number ? work_eta[i] = ηnm[i] : work_eta[i] = ηnm[i](t)
            end
            return work_eta
        end
        if displacement == "truncated" && !time_dependent_eta
            D_arrays = []
            for (i, mode) in enumerate(allmodes)
                push!(D_arrays, real.(displace(mode, ηlist(0)[i]).data))
            end
        end

        # iterate over ion-laser transitions
        for (ti, tr) in enumerate(ts)
            Δ, Ω = Δm[rn, m][ti], Ωm[rn, m][ti]
            Δ_2π = Δ / 2π
            typeof(Ω) <: Number && continue  # e.g. the laser doesn't shine on this ion
            locs = view(ion_idxs, ((ti-1)*ion_reps+1):(ti*ion_reps))
            for j in 1:prod(mode_dims)
                for i in nzrange(ld_array, j)
                    ri = rows[i]
                    ri < j && continue
                    cf = vals[i]
                    pflag = abs(Δ_2π + cf) > rwa_cutoff
                    nflag = abs(Δ_2π - cf) > rwa_cutoff
                    (pflag && nflag) && continue
                    rev_indxs = false
                    idxs = _inv_get_kron_indxs((rows[i], j), mode_dims)
                    for l in 1:L
                        (idxs[1][l] ≠ idxs[2][l] && typeof(ηnm[l]) <: Number) && @goto cl
                    end
                    s_ri = []
                    s_cri = []
                    for loc in locs
                        if !pflag
                            push!(
                                s_ri,
                                (Q * (ri - 1) + ion_rows[loc], Q * (j - 1) + ion_cols[loc])
                            )
                            if !nflag && ri ≠ j
                                push!(
                                    s_cri,
                                    (
                                        Q * (j - 1) + ion_rows[loc],
                                        Q * (ri - 1) + ion_cols[loc]
                                    )
                                )
                            end
                        elseif !nflag
                            push!(
                                s_ri,
                                (Q * (j - 1) + ion_rows[loc], Q * (ri - 1) + ion_cols[loc])
                            )
                            rev_indxs = true
                        end
                    end
                    if rev_indxs
                        idxs = reverse(idxs)
                    end
                    if length(s_cri) > 0
                        parity = sum(map(x -> isodd(abs(x[1] - x[2])), zip(idxs...)))
                        pushfirst!(s_cri, (-1 * isodd(parity), 0))
                    end

                    # push information to top-level lists/ construct time-dep function
                    if displacement == "truncated" && !time_dependent_eta
                        D = Tuple([D_arrays[i][idxs[1][i], idxs[2][i]] for i in 1:L])
                    elseif displacement == "analytic" && !time_dependent_eta
                        D = Tuple([
                            _Dnm_cnst_eta(ηlist(0)[i], idxs[1][i], idxs[2][i]) for i in 1:L
                        ])
                    elseif displacement == "truncated"
                        pflist = [_pf(mode_dims[i], idxs[1][i], idxs[2][i]) for i in 1:L]
                    end
                    row, col = s_ri[1]
                    if haskey(indxs_dict, s_ri[1])
                        # this will happen when multiple lasers address the same transition
                        functions[indxs_dict[row, col]] = let
                            a = functions[indxs_dict[row, col]]
                            if !time_dependent_eta
                                FunctionWrapper{Tuple{ComplexF64, ComplexF64}, Tuple{Float64}}(
                                    t ->
                                        a(t) .+ _D_cnst_eta(
                                            Ω(t),
                                            Δ,
                                            νlist,
                                            timescale,
                                            idxs,
                                            D,
                                            t,
                                            L
                                        )
                                )
                            elseif displacement == "analytic"
                                FunctionWrapper{Tuple{ComplexF64, ComplexF64}, Tuple{Float64}}(
                                    t ->
                                        a(t) .+ _D(
                                            Ω(t),
                                            Δ,
                                            ηlist(t),
                                            νlist,
                                            timescale,
                                            idxs,
                                            t,
                                            L
                                        )
                                )
                            elseif displacement == "truncated"
                                FunctionWrapper{Tuple{ComplexF64, ComplexF64}, Tuple{Float64}}(
                                    t ->
                                        a(t) .+ _Dtrunc(
                                            Ω(t),
                                            Δ,
                                            ηlist(t),
                                            νlist,
                                            rootlist,
                                            mode_dims,
                                            idxs,
                                            pflist,
                                            timescale,
                                            L,
                                            t
                                        )
                                )
                            end
                        end
                    else
                        if !time_dependent_eta
                            f = FunctionWrapper{
                                Tuple{ComplexF64, ComplexF64},
                                Tuple{Float64}
                            }(
                                t -> _D_cnst_eta(Ω(t), Δ, νlist, timescale, idxs, D, t, L)
                            )
                        elseif displacement == "analytic"
                            f = FunctionWrapper{
                                Tuple{ComplexF64, ComplexF64},
                                Tuple{Float64}
                            }(
                                t ->
                                    _D(Ω(t), Δ, ηlist(t), νlist, timescale, idxs, t, L)
                            )
                        elseif displacement == "truncated"
                            f = FunctionWrapper{
                                Tuple{ComplexF64, ComplexF64},
                                Tuple{Float64}
                            }(
                                t -> _Dtrunc(
                                    Ω(t),
                                    Δ,
                                    ηlist(t),
                                    νlist,
                                    rootlist,
                                    mode_dims,
                                    idxs,
                                    pflist,
                                    timescale,
                                    L,
                                    t
                                )
                            )
                        end
                        push!(functions, f)
                        push!(repeated_indices, s_ri)
                        push!(conj_repeated_indices, s_cri)
                        indxs_dict[row, col] = length(repeated_indices)
                    end
                    @label cl
                end
            end
        end
    end
    return functions, repeated_indices, conj_repeated_indices
end

# δν(t) × aᵀa terms for Hamiltonian. This function returns an array of functions
# δν_functions = [2π×ν.δν(t)×timescale for ν in modes]. It also returns an array of arrays
# of arrays of indices, δν_indices, such that δν_indices[i][j] lists all diagonal elements
# of the full 2D system matrix upon which have been mapped the jth diagonal element of the
# ith mode.
function _setup_δν_hamiltonian(chamber, timescale)
    N = length(ions(chamber))
    allmodes = modes(chamber)
    δν_indices = Vector{Vector{Vector{Int64}}}(undef, 0)
    δν_functions = FunctionWrapper[]
    τ = timescale
    for l in eachindex(allmodes)
        δν = frequency_fluctuation(allmodes[l])
        mode = allmodes[l]
        (mode._cnst_δν && δν(0) == 0) && continue
        push!(
            δν_functions,
            FunctionWrapper{Float64, Tuple{Float64}}(t -> @fastmath 2π * δν(t) * τ)
        )
        δν_indices_l = Vector{Vector{Int64}}(undef, 0)
        mode_op = number(mode)
        A = embed(basis(chamber), [N + l], [mode_op]).data
        mode_dim = mode.shape[1]
        for i in 1:(mode_dim-1)
            indices = [x[1] for x in getfield.(findall(x -> x .== complex(i, 0), A), :I)]
            push!(δν_indices_l, indices)
        end
        push!(δν_indices, δν_indices_l)
    end
    return δν_indices, δν_functions
end

# Hamiltonian terms for global Bfield fluctuations encoded in T.δB. If T.δB=0, this function
# returns a collection of empty arrays. Otherwise it iterates over the selected levels for
# each ion and creates an array (global_B_scales), which encodes the magnetic field
# susceptibility of each level. It also returns an array of indices (global_B_indices), that
# keeps track of the indices of the full 2D array that represents the tensored system,
# corresponding to the energy of each level. Finally it returns a function (bfunc) encoding
# the time-dependence of δB. When the system is integrated, the Hamiltonian terms will be
# updated at each time step by by bfunc(t) times the individual susceptibilities.
function _setup_global_B_hamiltonian(chamber, timescale)
    all_ions = ions(chamber)
    global_B_indices = Vector{Vector{Int64}}(undef, 0)
    global_B_scales = Vector{Float64}(undef, 0)
    δB = magnitude_fluctuation(magneticfield(chamber))
    τ = timescale
    bfunc = FunctionWrapper{Float64, Tuple{Float64}}(t -> 2π * δB(t * τ))
    if chamber._cnst_δB && δB(0) == 0
        return global_B_indices, global_B_scales, bfunc
    end
    for n in eachindex(all_ions)
        for sublevel in sublevels(all_ions[n])
            ion_op = sigma(all_ions[n], sublevel)
            A = embed(basis(chamber), [n], [ion_op]).data
            indices = [x[1] for x in getfield.(findall(x -> x .== complex(1, 0), A), :I)]
            push!(global_B_indices, indices)
            # zeemanshift(ions[n], sublevel, 1]) is the Zeeman shift of
            # sublevel in units of δB.
            push!(global_B_scales, τ * zeemanshift(all_ions[n], sublevel, 1))
        end
    end
    return global_B_indices, global_B_scales, bfunc
end

# This mostly just strings together the results from _setup_global_B_hamiltonian and
# _setup_δν_hamiltonian for use in the hamiltonian function. The one additional task
# performed is the creation of an array of indices (all_unique_indices), which keeps track
# of all the diagonal indices affected by δν and/or δB, which is useful in hamiltonian().
function _setup_fluctuation_hamiltonian(chamber, timescale)
    gbi, gbs, bfunc = _setup_global_B_hamiltonian(chamber, timescale)
    δνi, δνfuncs = _setup_δν_hamiltonian(chamber, timescale)
    all_unique_indices = convert(Vector{Int64}, _flattenall(unique([gbi; δνi])))
    return all_unique_indices, gbi, gbs, bfunc, δνi, δνfuncs
end

#############################################################################################
# internal functions
#############################################################################################

# https://gist.github.com/ivirshup/e9148f01663278ca4972d8a2d9715f72
function _flattenall(a::AbstractArray)
    while any(x -> typeof(x) <: AbstractArray, a)
        a = collect(Iterators.flatten(a))
    end
    return a
end

# A 3D array of Lamb-Dicke parameters for each combination of ion, laser and mode. Modes are
# populated in reverse order.
function _ηmatrix(T)
    all_ions = ions(T)
    vms = modes(T)
    all_lasers = lasers(T)
    (N, M, L) = map(x -> length(x), [all_ions, all_lasers, vms])
    ηnml = Array{Any}(undef, N, M, L)
    for n in 1:N, m in 1:M, l in 1:L
        δν = frequency_fluctuation(vms[l])
        ν = frequency(vms[l])
        eta = lambdicke(vms[l], all_ions[n], all_lasers[m], scaled=true)
        if eta == 0
            ηnml[n, m, L-l+1] = 0
        else
            ηnml[n, m, L-l+1] =
                FunctionWrapper{Float64, Tuple{Float64}}(t -> eta / √(ν + δν(t)))
        end
    end
    return ηnml
end

# Returns an array of vectors. The rows and columns of the array refer to ions and lasers,
# respectively. For each row/column we have a vector of detunings from the laser frequency
# for each ion transition. We need to separate this calculation from _Ωmatrix to implement
# RWA easily.
function _Δmatrix(chamber, timescale)
    all_ions = ions(chamber)
    all_lasers = lasers(chamber)
    (N, M) = length(all_ions), length(all_lasers)
    Δnmkj = Array{Vector}(undef, N, M)
    for n in 1:N, m in 1:M
        Btot = bfield(chamber, all_ions[n])
        v = Vector{Float64}(undef, 0)
        for transition in subleveltransitions(all_ions[n])
            ωa = transitionfrequency(all_ions[n], transition, B=Btot)
            push!(
                v,
                2π *
                timescale *
                ((c / wavelength(all_lasers[m])) + detuning(all_lasers[m]) - ωa)
            )
        end
        Δnmkj[n, m] = v
    end
    return Δnmkj
end

# Returns an array of vectors. the rows and columns of the array refer to ions and lasers,
# respectively. For each row/column we have a vector of coupling strengths between the laser
# and all allowed electronic ion transitions.
function _Ωmatrix(chamber, timescale)
    all_ions = ions(chamber)
    all_lasers = lasers(chamber)
    (N, M) = length(all_ions), length(all_lasers)
    Ωnmkj = Array{Vector}(undef, N, M)
    for n in 1:N, m in 1:M
        I = intensity(all_lasers[m])
        ϕ = phase(all_lasers[m])
        transitions = subleveltransitions(all_ions[n])
        s_indx = findall(x -> x[1] == n, pointing(all_lasers[m]))
        if length(s_indx) == 0
            Ωnmkj[n, m] = [0 for _ in 1:length(transitions)]
            continue
        else
            s = pointing(all_lasers[m])[s_indx[1]][2]
        end
        v = []
        for t in transitions
            Ω0 =
                2π *
                timescale *
                s *
                matrixelement(
                    all_ions[n],
                    t,
                    1.0,
                    polarization(all_lasers[m]),
                    wavevector(all_lasers[m]),
                    direction(magneticfield(chamber))
                ) / 2.0
            if Ω0 == 0
                push!(v, 0)
            else
                push!(
                    v,
                    FunctionWrapper{ComplexF64, Tuple{Float64}}(
                        t -> Ω0 * √I(t) * exp(-im * ϕ(t))
                    )
                )
            end
        end
        Ωnmkj[n, m] = v
    end
    return Ωnmkj
end

# Returns a tuple correpsonding to: [σ₊(t)]_ij ⋅ [D(ξ(t))]_ij, [σ₊(t)]_ji ⋅ [D(ξ(t))]_ji.
# [D(ξ(t))]_ij is calculated assuming an infinite dimensional Hilbert space for the HO.
function _D(Ω, Δ, η, ν, timescale, n, t, L)
    d = complex(1, 0)
    for i in 1:L
        d *= _Dnm(1im * η[i] * exp(im * 2π * ν[i] * timescale * t), n[1][i], n[2][i])
    end
    g = Ω * exp(-1im * t * Δ)
    return g * d, g * conj(d)
end

# Returns a tuple correpsonding to: [σ₊(t)]_ij ⋅ [D(ξ(t))]_ij, [σ₊(t)]_ji ⋅ [D(ξ(t))]_ji.
# [D(ξ(t))]_ij is calculated assuming an infinite dimensional Hilbert space for the HO.
# As opposed to _D, in this case, we assume η(t) = η₀, which allows us to precompute _Dnm.
# This precomputation is performed externally to the function and fed in as the argument `D`.
function _D_cnst_eta(Ω, Δ, ν, timescale, n, D, t, L)
    d = complex(1, 0)
    for i in 1:L
        d *= D[i] * exp(1im * (n[1][i] - n[2][i]) * (2π * ν[i] * timescale * t + π / 2))
    end
    g = Ω * exp(-1im * t * Δ)
    return g * d, g * conj(d)
end

# Consider: T = X₁ ⊗ X₂ ⊗ ... ⊗ X_n (Xᵢ ∈ ℝ{dims[i]×dims[i]}), and indices:
# indxs[1], indxs[2], ..., indsx[N] = (i1, j1), (i2, j2), ..., (iN, jN).
# This function returns (k, l) such that: T[k, l] = X₁[i1, j1] * X₂[i2, j2] *...* X_N[iN, jN]
function _get_kron_indxs(indxs::Vector{Tuple{Int64, Int64}}, dims::Vector{Int64})
    L = length(indxs)
    rowcol = Int64[0, 0]
    @simd for i in 0:(L-1)
        if i == 0
            @inbounds rowcol .+= indxs[L-i]
        else
            @inbounds rowcol .+= (indxs[L-i] .- 1) .* prod(view(dims, 1:i))
        end
    end
    return rowcol
end

# The inverse of _get_kron_indxs. If T = X₁ ⊗ X₂ ⊗ X₃ and X₁, X₂, X₃ are M×M, N×N and L×L
# dimension matrices, then we should input dims=(M, N, L).
function _inv_get_kron_indxs(indxs, dims)
    row, col = indxs
    N = length(dims)
    ret_rows = Array{Int64}(undef, N)
    ret_cols = Array{Int64}(undef, N)
    for i in 1:N
        tensor_N = prod(dims[i:N])
        M = tensor_N ÷ dims[i]
        rowflag = false
        colflag = false
        for j in 1:dims[i]
            jM = j * M
            if !rowflag && row <= jM
                @inbounds ret_rows[i] = j
                row -= jM - M
                rowflag = true
            end
            if !colflag && col <= jM
                @inbounds ret_cols[i] = j
                col -= jM - M
                colflag = true
            end
            rowflag && colflag && break
        end
    end
    return Tuple(ret_rows), Tuple(ret_cols)
end

# similar to _Dnm, but meant to be used when η is assumed constant in ξ=iηe^(i2πνt)
function _Dnm_cnst_eta(ξ::Number, n::Int, m::Int)
    n < m && return _Dnm_cnst_eta(ξ, m, n) * (-1)^isodd(abs(n - m))
    n -= 1
    m -= 1
    s = 1.0
    for i in (m+1):n
        s *= i
    end
    ret = sqrt(1 / s) * ξ^(n - m) * exp(-abs2(ξ) / 2.0) * _alaguerre(abs2(ξ), m, n - m)
    isnan(ret) && return 1.0 * (n == m)
    return ret
end

# If lamb_dicke_order is <: Int, this constructs a constant vector with this value of length
# L (i.e. same lamb_dicke_order for all modes). Otherwise lamb_dicke_order is reversed and
# returned.
function _check_lamb_dicke_order(lamb_dicke_order, L)
    if typeof(lamb_dicke_order) <: Int
        return [lamb_dicke_order for _ in 1:L]
    else
        @assert(
            length(lamb_dicke_order) == L,
            "if typeof(lamb_dicke_order)<:Vector, then length of lamb_dicke_order must ",
            "equal number of modes"
        )
        reverse(lamb_dicke_order)
    end
end

function _ld_array(mode_dims, lamb_dicke_order, νlist, timescale)
    a = [spzeros(Float16, d, d) for d in mode_dims]
    @inbounds for (i, d) in enumerate(mode_dims)
        for k in 1:d, l in 1:k
            if k - l <= lamb_dicke_order[i]
                val = (l - k) * νlist[i] * timescale
                a[i][k, l] = exp(val)
                l ≠ k && @inbounds a[i][l, k] = exp(-val)
            end
        end
    end
    length(a) == 1 ? ld_array = a[1] : ld_array = kron(a...)
    return ld_array, rowvals(ld_array), log.(nonzeros(ld_array))
end
