using SparseArrays: rowvals, nzrange, nonzeros
using FunctionWrappers: FunctionWrapper
using PolynomialRoots: roots
using QuantumOptics: SparseOperator, embed

export hamiltonian


"""
    hamiltonian(
            T::Trap; timescale::Real=1e-6, lamb_dicke_order::Union{Vector{Int},Int}=1, 
            rwa_cutoff::Real=Inf, displacement="truncated", time_dependent_eta=false
        )      
Constructs the Hamiltonian for `T` as a function of time. Return type is a function 
`h(t::Real, ψ)` that, itself, returns a `QuantumOptics.SparseOperator`.

#### args
* `timescale`: e.g. a value of 1e-6 will take time to be in ``\\mu s``
* `lamb_dicke_order`: Only consider terms that change the phonon number by up to this value.
    If this is an `Int`, then the cutoff is applied to all modes. If this is a `Vector{Int}`,
    then `lamb_dicke_order[i]` is applied to the iᵗʰ mode, according to the order in 
    `T.basis`.
    Note: this isn't quite the same thing as the Lamb-Dicke approximation since setting
    `lamb_dicke_order=1` will retain, for example, terms proportional to ``a^\\dagger a ``.
* `rwa_cutoff`: drop terms in the Hamiltonian that oscillate faster than this cutoff. **Note:
    if not using an RWA set to `Inf` (rather than a large number) for faster performance.**
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
        T::Trap; timescale::Real=1e-6, lamb_dicke_order::Union{Vector{Int},Int}=1, 
        rwa_cutoff::Real=Inf, displacement::String="truncated", time_dependent_eta::Bool=false
    ) 
    hamiltonian(
            T, T.configuration, timescale, lamb_dicke_order, rwa_cutoff, displacement,
            time_dependent_eta
        ) 
end


#############################################################################################
# Hamiltonian for a linear configuration of ions
#############################################################################################

# At each time step, this function updates in-place the 2D array describing the full system
# Hamiltonian. 
function hamiltonian(
        T::Trap, configuration::LinearChain, timescale::Real, 
        lamb_dicke_order::Union{Vector{Int},Int}, rwa_cutoff::Real, displacement::String,
        time_dependent_eta::Bool
    )
    b, indxs, cindxs = _setup_base_hamiltonian(
            T, timescale, lamb_dicke_order, rwa_cutoff, displacement, time_dependent_eta
        )
    aui, gbi, gbs, bfunc, δνi, δνfuncs = _setup_fluctuation_hamiltonian(T, timescale)
    S = SparseOperator(get_basis(T))
    function f(t, ψ)  # a two argument function is required in the QuantumOptics solvers
        @inbounds begin
            @simd for i in 1:length(indxs)
                bt_i, conj_bt_i = b[i](t)::Tuple{ComplexF64,ComplexF64}
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
                    S.data[indx, indx] = complex(0)
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

The _setup_hamiltonian function extracts this user-defined information from a <:Trap struct 
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
        T, timescale, lamb_dicke_order, rwa_cutoff, displacement, time_dependent_eta
    )
    modes = reverse(get_vibrational_modes(T.configuration))
    L = length(modes)
    ηm, Δm, Ωm = _ηmatrix(T), _Δmatrix(T, timescale), _Ωmatrix(T, timescale)
    ions = T.configuration.ions
    indxs_dict = Dict()
    repeated_indices = Vector{Vector{Tuple{Int64,Int64}}}(undef, 0)
    conj_repeated_indices = Vector{Vector{Tuple{Int64,Int64}}}(undef, 0)
    functions = FunctionWrapper[]
    νlist = [mode.ν for mode in modes]
    mode_dims = [mode.N+1 for mode in modes]
    N = prod(mode_dims)
    lamb_dicke_order = _check_lamb_dicke_order(lamb_dicke_order, L)
    ld_arrays, rootlist = _construct_ld_arrays(
        mode_dims, lamb_dicke_order, νlist, timescale, displacement, time_dependent_eta
    )
    νarray = zeros(Float64, N, N)
    indx_array = zeros(ComplexF64, N, N)
    work_eta = zeros(Float64, L)

    # iterate over ions, lasers and ion-laser transitions
    for n in eachindex(ions), m in eachindex(T.lasers), (ti, tr) in enumerate(_transitions(ions[n]))
        ηnm = view(ηm, n, m, :)
        function ηlist(t)
            for i in 1:L
                work_eta[i] = ηnm[i](t)
            end
            work_eta
        end
        Δ, Ω = Δm[n, m][ti], Ωm[n, m][ti]
        if typeof(Ω) <: Number
            # e.g. the laser doesn't shine on this ion
            continue 
        end

        # Construct a complex array where each element has a value z, with real(z) equal to 
        # the corresponding row and imag(z) equal to the corresponding column. We will 
        # construct a kronecker product with this array taking the place of a deisplacement 
        # operator and then use z to determine where the elements have been mapped. Also, 
        # performs Lamb-Dicke approx. + RWA by constructing an array of indices with nonzero 
        # values only when matrix element satisfies both. 
        νarray .*= 0.; indx_array .*= 0.
        if length(ld_arrays) > 1
            νarray .= kron(
                [typeof(ηnm[l])<:Number ? one(ld_arrays[l]) : ld_arrays[l] for l in 1:L]...
            )
        else
            νarray .= ld_arrays[1]
        end
        for i in 1:N, j in 1:N
            if νarray[i, j] == 0
                continue
            elseif abs((Δ/(2π)) + log(νarray[i, j])) < rwa_cutoff * timescale
                indx_array[i, j] = complex(i, j)
            end
        end

        # construct the tensor product 𝐼 ⊗...⊗ σ₊ ⊗...⊗ 𝐼 ⊗ indx_array
        ion_op = sigma(ions[n], tr[2], tr[1])
        mode_op = SparseOperator(⊗(reverse(modes)...), indx_array)
        I = length(ions)
        A = embed(get_basis(T), [n, collect(I+1:I+L)], [ion_op, mode_op]).data

        if displacement == "truncated" && !time_dependent_eta
            D_arrays = []
            for (i, mode) in enumerate(modes)
                push!(D_arrays, real.(displace(mode, ηlist(0)[i]).data))
            end
        end

        # See where subspace operators have been mapped after embedding
        for i in 1:N, j in (isinf(rwa_cutoff) ? collect(1:i) : 1:N)
            sub_idxs = _inv_get_kron_indxs([i, j], mode_dims)
            # find all locations of i + im*j 
            s_ri = sort(getfield.(findall(x->x.==complex(i, j), A), :I), by=x->x[2])
            if length(s_ri) == 0  
                # this index doesn't exist due to Lamd-Dicke approx and/or RWA
                continue
            end
            # if not RWA, find all locations of j, i
            if isinf(rwa_cutoff) && i != j
                s_cri = sort(getfield.(findall(x->x.==complex(j, i), A), :I), by=x->x[2])
                parity = sum(map(x->isodd(abs(x[1]-x[2])), zip(sub_idxs...)))
                pushfirst!(s_cri, (-1 * isodd(parity), 0))
            else
                s_cri = []
            end

            # push information to top-level lists/ construct time-dep function
            if displacement == "truncated" && !time_dependent_eta
                D = [D_arrays[i][sub_idxs[1][i], sub_idxs[2][i]] for i in 1:L]
            elseif displacement == "analytic" && !time_dependent_eta
                D = [_Dnm_cnst_eta(ηlist(0)[i], sub_idxs[1][i], sub_idxs[2][i]) for i in 1:L]
            elseif displacement == "truncated"
                pflist = [_pf(mode_dims[i], sub_idxs[1][i], sub_idxs[2][i]) for i in 1:L]
            end
            row, col = s_ri[1]
            if haskey(indxs_dict, s_ri[1]) 
                # this will happen when multiple lasers address the same transition
                functions[indxs_dict[row, col]] = 
                    let 
                        a = functions[indxs_dict[row, col]]
                        if !time_dependent_eta
                            FunctionWrapper{Tuple{ComplexF64,ComplexF64},Tuple{Float64}}(
                                t -> a(t) .+ _D_cnst_eta(
                                    Ω(t), Δ, νlist, timescale, sub_idxs, D, t, L))
                        elseif displacement == "analytic"
                            FunctionWrapper{Tuple{ComplexF64,ComplexF64},Tuple{Float64}}(
                                t -> a(t) .+ _D(
                                    Ω(t), Δ, ηlist(t), νlist, timescale, sub_idxs, t, L))
                        elseif displacement == "truncated"
                            FunctionWrapper{Tuple{ComplexF64,ComplexF64},Tuple{Float64}}(
                                t -> a(t) .+ _Dtrunc( 
                                    Ω(t), Δ, ηlist(t), νlist, rootlist, mode_dims, sub_idxs,
                                    pflist, timescale, L, t))
                        end
                    end
            else
                if !time_dependent_eta
                    f = FunctionWrapper{Tuple{ComplexF64,ComplexF64},Tuple{Float64}}(
                            t -> _D_cnst_eta(Ω(t), Δ, νlist, timescale, sub_idxs, D, t, L))
                elseif displacement == "analytic"
                    f = FunctionWrapper{Tuple{ComplexF64,ComplexF64},Tuple{Float64}}(
                            t -> _D(Ω(t), Δ, ηlist(t), νlist, timescale, sub_idxs, t, L))
                elseif displacement == "truncated"
                    f = FunctionWrapper{Tuple{ComplexF64,ComplexF64},Tuple{Float64}}(
                                t -> _Dtrunc( 
                                    Ω(t), Δ, ηlist(t), νlist, rootlist, mode_dims, sub_idxs,
                                    pflist, timescale, L, t))
                end
                push!(functions, f)
                push!(repeated_indices, s_ri)
                push!(conj_repeated_indices, s_cri)
                indxs_dict[row, col] = length(repeated_indices)
            end
        end
    end
    functions, repeated_indices, conj_repeated_indices
end

# δν(t) × aᵀa terms for Hamiltonian. This function returns an array of functions
# δν_functions = [2π×ν.δν(t)×timescale for ν in modes]. It also returns an array of arrays
# of arrays of indices, δν_indices, such that δν_indices[i][j] lists all diagonal elements
# of the full 2D system matrix upon which have been mapped the jth diagonal element of the 
# ith mode. 
function _setup_δν_hamiltonian(T, timescale)
    N = length(T.configuration.ions)
    modes = get_vibrational_modes(T.configuration)
    δν_indices = Vector{Vector{Vector{Int64}}}(undef, 0)
    δν_functions = FunctionWrapper[]
    τ = timescale
    for l in eachindex(modes)
        δν = modes[l].δν
        mode = modes[l]
        if mode._cnst_δν && δν(0) == 0
            continue
        end
        push!(δν_functions, FunctionWrapper{Float64,Tuple{Float64}}(
                    t -> @fastmath 2π * δν(t) * τ
                ))
        δν_indices_l = Vector{Vector{Int64}}(undef, 0)
        mode_op = number(mode)
        A = embed(get_basis(T), [N+l], [mode_op]).data
        mode_dim = mode.shape[1]
        for i in 1:mode_dim-1
            indices = [x[1] for x in getfield.(findall(x->x.==complex(i, 0), A), :I)]
            push!(δν_indices_l, indices)
        end
        push!(δν_indices, δν_indices_l)
    end
    δν_indices, δν_functions
end

# Hamiltonian terms for global Bfield fluctuations encoded in T.δB. If T.δB=0, this function
# returns a collection of empty arrays. Otherwise it iterates over the selected levels for
# each ion and creates an array (global_B_scales), which encodes the magnetic field
# susceptibility of each level. It also returns an array of indices (global_B_indices), that
# keeps track of the indices of the full 2D array that represents the tensored system, 
# corresponding to the energy of each level. Finally it returns a function (bfunc) encoding 
# the time-dependence of δB. When the system is integrated, the Hamiltonian terms will be 
# updated at each time step by by bfunc(t) times the individual susceptibilities.
function _setup_global_B_hamiltonian(T, timescale)
    ions = T.configuration.ions
    global_B_indices = Vector{Vector{Int64}}(undef, 0)
    global_B_scales = Vector{Float64}(undef, 0)
    δB = T.δB
    τ = timescale
    bfunc = FunctionWrapper{Float64,Tuple{Float64}}(t -> 2π * δB(t) * τ)
    if T._cnst_δB && δB(0) == 0
        return global_B_indices, global_B_scales, bfunc
    end
    for n in eachindex(ions)
        transitions_list = [] 
        for (ti, tr) in enumerate(_transitions(ions[n]))
            for level in tr
                if level in transitions_list
                    continue
                else
                    push!(transitions_list, level)
                end
                ion_op = sigma(ions[n], level)
                A = embed(get_basis(T), [n], [ion_op]).data
                indices = [x[1] for x in getfield.(findall(x->x.==complex(1, 0), A), :I)]
                push!(global_B_indices, indices)
                # zeeman_shift(1, ions[n].level_structure[level]) is the Zeeman shift of
                # ions[n].level_structure[level] in units of δB.
                push!(global_B_scales, τ * zeeman_shift(1, ions[n].level_structure[level]))
            end
        end
    end
    global_B_indices, global_B_scales, bfunc
end

# This mostly just strings together the results from _setup_global_B_hamiltonian and 
# _setup_δν_hamiltonian for use in the hamiltonian function. The one additional task 
# performed is the creation of an array of indices (all_unique_indices), which keeps track
# of all the diagonal indices affected by δν and/or δB, which is useful in hamiltonian().
function _setup_fluctuation_hamiltonian(T, timescale)
    gbi, gbs, bfunc = _setup_global_B_hamiltonian(T, timescale)
    δνi, δνfuncs = _setup_δν_hamiltonian(T, timescale)
    all_unique_indices = convert(Vector{Int64}, _flattenall(unique([gbi; δνi])))
    all_unique_indices, gbi, gbs, bfunc, δνi, δνfuncs
end


#############################################################################################
# internal functions
#############################################################################################

# https://gist.github.com/ivirshup/e9148f01663278ca4972d8a2d9715f72
function _flattenall(a::AbstractArray)
    while any(x -> typeof(x)<:AbstractArray, a)
        a = collect(Iterators.flatten(a))
    end
    a
end

# A 3D array of Lamb-Dicke parameters for each combination of ion, laser and mode. Modes are
# populated in reverse order.
function _ηmatrix(T)
    ions = T.configuration.ions
    vms = get_vibrational_modes(T.configuration)
    lasers = T.lasers
    (N, M, L) = map(x -> length(x), [ions, lasers, vms])
    ηnml = Array{Any}(undef, N, M, L)
    for n in 1:N, m in 1:M, l in 1:L
        δν = vms[l].δν
        ν = vms[l].ν
        eta = get_η(vms[l], lasers[m], ions[n], scaled=true)
        if eta == 0
            ηnml[n,m,L-l+1] = 0
        else
            ηnml[n,m,L-l+1] = FunctionWrapper{Float64,Tuple{Float64}}(t->eta / √(ν + δν(t)))
        end
    end
    ηnml
end

# Returns an array of vectors. The rows and columns of the array refer to ions and lasers,
# respectively. For each row/column we have a vector of detunings from the laser frequency 
# for each ion transition. We need to separate this calculation from  _Ωmatrix to implement
# RWA easily.  
function _Δmatrix(T, timescale)
    ions = T.configuration.ions
    lasers = T.lasers
    (N, M) = length(ions), length(lasers)
    B = T.B
    ∇B = T.∇B
    Δnmkj = Array{Vector}(undef, N, M)
    for n in 1:N, m in 1:M
        transitions = _transitions(ions[n])
        Btot = B + ∇B * ions[n].position
        v = Vector{Float64}(undef, 0)
        for t in transitions
            L1 = ions[n].selected_level_structure[t[1]]
            L2 = ions[n].selected_level_structure[t[2]]
            stark_shift = ions[n].stark_shift[t[1]] - ions[n].stark_shift[t[2]]
            ωa = (abs(L1.E  + zeeman_shift(Btot, L1) - (L2.E + zeeman_shift(Btot, L2)))
                  + stark_shift)
            push!(v, 2π * timescale * ((c / lasers[m].λ) + lasers[m].Δ - ωa))
        end
        Δnmkj[n, m] = v
    end
    Δnmkj
end

# Returns an array of vectors. the rows and columns of the array refer to ions and lasers,
# respectively. For each row/column we have a vector of coupling strengths between the laser
# and all allowed electronic ion transitions.
function _Ωmatrix(T, timescale)
    ions = T.configuration.ions
    lasers = T.lasers
    (N, M) = length(ions), length(lasers)
    Ωnmkj = Array{Vector}(undef, N, M)
    for n in 1:N, m in 1:M
        E = lasers[m].E
        phase = lasers[m].ϕ
        (γ, ϕ) = map(x -> rad2deg(acos(ndot(T.Bhat, x))), [lasers[m].ϵ, lasers[m].k])
        transitions = _transitions(ions[n])
        s_indx = findall(x -> x[1] == n, lasers[m].pointing) 
        if length(s_indx) == 0
            Ωnmkj[n, m] = [0 for _ in 1:length(transitions)]
            continue
        else
            s = lasers[m].pointing[s_indx[1]][2]
        end
        v = []
        for t in transitions
            Ω0 = 2π * timescale * s * ions[n].selected_matrix_elements[t](1.0, γ, ϕ) / 2.0
            if Ω0 == 0
                push!(v, 0)
            else
                push!(v, FunctionWrapper{ComplexF64,Tuple{Float64}}(
                        t -> Ω0 * E(t) * exp(-im * phase(t)))
                    )
            end
        end
        Ωnmkj[n, m] = v
    end
    Ωnmkj
end

_transitions(ion) = collect(keys(ion.selected_matrix_elements))

# Returns a tuple correpsonding to: [σ₊(t)]_ij ⋅ [D(ξ(t))]_ij, [σ₊(t)]_ji ⋅ [D(ξ(t))]_ji.
# [D(ξ(t))]_ij is calculated assuming an infinite dimensional Hilbert space for the HO.
function _D(Ω, Δ, η, ν, timescale, n, t, L)
    d = complex(1, 0)
    @simd for i in 1:L
        d *= _Dnm(1im * η[i] * exp(im * 2π * ν[i] * timescale * t), n[1][i], n[2][i])
    end
    g = Ω * exp(-1im * t * Δ)
    g * d, g * conj(d)
end

# Returns a tuple correpsonding to: [σ₊(t)]_ij ⋅ [D(ξ(t))]_ij, [σ₊(t)]_ji ⋅ [D(ξ(t))]_ji.
# [D(ξ(t))]_ij is calculated assuming an infinite dimensional Hilbert space for the HO.
# As opposed to _D, in this case, we assume η(t) = η₀, which allows us to precompute _Dnm. 
# This precomputation is performed externally to the function and fed in as the argument `D`.
function _D_cnst_eta(Ω, Δ, ν, timescale, n, D, t, L)
    d = complex(1, 0)
    @simd for i in 1:L
        d *= D[i] * exp(1im * (n[1][i] - n[2][i]) * (2π * ν[i] * timescale * t + π/2))
    end
    g = Ω * exp(-1im * t * Δ)
    g * d, g * conj(d)
end

# Consider: T = X₁ ⊗ X₂ ⊗ ... ⊗ X_n (Xᵢ ∈ ℝ{dims[i]×dims[i]}), and indices: 
# indxs[1], indxs[2], ..., indsx[N] = (i1, j1), (i2, j2), ..., (iN, jN). 
# This function returns (k, l) such that: T[k, l] = X₁[i1, j1] * X₂[i2, j2] *...* X_N[iN, jN]
function _get_kron_indxs(indxs, rdims)
    dims = reverse(rdims)
    row, col = 0, 0
    for (i, indx) in enumerate(reverse(indxs))
        if i == 1
            row += indx[1]
            col += indx[2]
        else
            row += (indx[1] - 1) * prod(dims[1:i-1])
            col += (indx[2] - 1) * prod(dims[1:i-1])
        end
    end
    row, col
end

# The inverse of _get_kron_indxs. If T = X₁ ⊗ X₂ ⊗ X₃ and X₁, X₂, X₃ are M×M, N×N and L×L
# dimension matrices, then we should input dims=(M, N, L). 
function _inv_get_kron_indxs(indxs, dims)
    row, col = indxs
    N = length(dims)
    ret_rows, ret_cols = Int64[], Int64[]
    for i in 1:N
        tensor_N = prod(dims[i:N])
        M = tensor_N ÷ dims[i]
        rowflag, colflag = false, false
        for j in 1:dims[i]
            if !rowflag && row <= j * M
                push!(ret_rows, j)
                row -= (j - 1) * M
                rowflag = true
            end
            if !colflag && col <= j * M 
                push!(ret_cols, j)
                col -= (j - 1) * M
                colflag = true
            end
            if rowflag && colflag
                break
            end
        end
    end
    tuple(ret_rows...), tuple(ret_cols...)
end

# similar to _Dnm, but meant to be used when η is assumed constant in ξ=iηe^(i2πνt)
function _Dnm_cnst_eta(ξ::Number, n::Int, m::Int)
    if n < m 
        return _Dnm_cnst_eta(ξ, m, n) * (-1)^isodd(abs(n-m))
    end
    n -= 1; m -= 1
    s = 1.0
    for i in m+1:n
        s *= i
    end
    ret = sqrt(1 / s) * ξ^(n-m) * exp(-abs2(ξ) / 2.0) * _alaguerre(abs2(ξ), m, n-m)
    if isnan(ret)
        return 1.0 * (n==m)
    end
    ret
end

# If lamb_dicke_order is <: Int, this constructs a constant vector with this value of length
# L. Otherwise lamb_dicke_order is reversed and returned.
function _check_lamb_dicke_order(lamb_dicke_order, L)
    if typeof(lamb_dicke_order) <: Int
        return [lamb_dicke_order for _ in 1:L]
    else
        @assert length(lamb_dicke_order) == L string(
            "if typeof(lamb_dicke_order)<:Vector, then length of lamb_dicke_order must ",
            "equal number of modes"
        )
        return reverse(lamb_dicke_order)
    end
end

function _construct_ld_arrays(
            mode_dims, lamb_dicke_order, νlist, timescale, displacement, time_dependent_eta
        )
    ld_arrays = []; rootlist = []
    for (i, dim) in enumerate(mode_dims)
        array = zeros(Float64, dim, dim)
        for k in 1:dim, l in 1:dim
            if abs(k - l) <= lamb_dicke_order[i]
                array[k, l] = exp((l - k) * νlist[i] * timescale)
            end
        end
        push!(ld_arrays, array)
        if displacement == "truncated" && time_dependent_eta
            push!(rootlist, real.(roots(_He(dim), polish=true)))
        end
    end
    ld_arrays, rootlist
end