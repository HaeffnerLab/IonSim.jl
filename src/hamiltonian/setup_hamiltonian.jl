
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
        (mode._cnst_δν && δν(0) == 0) && continue
        push!(
            δν_functions,
            FunctionWrapper{Float64, Tuple{Float64}}(t -> @fastmath 2π * δν(t) * τ)
        )
        δν_indices_l = Vector{Vector{Int64}}(undef, 0)
        mode_op = number(mode)
        A = embed(get_basis(T), [N + l], [mode_op]).data
        mode_dim = mode.shape[1]
        for i in 1:(mode_dim - 1)
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
function _setup_global_B_hamiltonian(T, timescale)
    ions = T.configuration.ions
    global_B_indices = Vector{Vector{Int64}}(undef, 0)
    global_B_scales = Vector{Float64}(undef, 0)
    δB = T.δB
    τ = timescale
    bfunc = FunctionWrapper{Float64, Tuple{Float64}}(t -> 2π * δB(t * τ))
    if T._cnst_δB && δB(0) == 0
        return global_B_indices, global_B_scales, bfunc
    end
    for n in eachindex(ions)
        for sublevel in sublevels(ions[n])
            ion_op = sigma(ions[n], sublevel)
            A = embed(get_basis(T), [n], [ion_op]).data
            indices = [x[1] for x in getfield.(findall(x -> x .== complex(1, 0), A), :I)]
            push!(global_B_indices, indices)
            # zeeman_shift(ions[n], sublevel, 1]) is the Zeeman shift of
            # sublevel in units of δB.
            push!(global_B_scales, τ * zeeman_shift(ions[n], sublevel, 1))
        end
    end
    return global_B_indices, global_B_scales, bfunc
end

# This mostly just strings together the results from _setup_global_B_hamiltonian and
# _setup_δν_hamiltonian for use in the hamiltonian function. The one additional task
# performed is the creation of an array of indices (all_unique_indices), which keeps track
# of all the diagonal indices affected by δν and/or δB, which is useful in hamiltonian().
function _setup_fluctuation_hamiltonian(T, timescale)
    gbi, gbs, bfunc = _setup_global_B_hamiltonian(T, timescale)
    δνi, δνfuncs = _setup_δν_hamiltonian(T, timescale)
    all_unique_indices = convert(Vector{Int64}, _flattenall(unique([gbi; δνi])))
    return all_unique_indices, gbi, gbs, bfunc, δνi, δνfuncs
end


"""
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
"""
function _setup_base_hamiltonian(
    T,
    timescale,
    lamb_dicke_order,
    rwa_cutoff,
    displacement,
    time_dependent_eta
)
    rwa_cutoff *= timescale
    modes = reverse(get_vibrational_modes(T.configuration))
    L = length(modes)
    νlist = Tuple([mode.ν for mode in modes])
    mode_dims = [mode.N + 1 for mode in modes]

    ions = reverse(T.configuration.ions)
    Q = prod([ion.shape[1] for ion in ions])
    ion_arrays = [spdiagm(0 => [true for _ in 1:ion.shape[1]]) for ion in ions]

    ηm, Δm, Ωm = _ηmatrix(T), _Δmatrix(T, timescale), _Ωmatrix(T, timescale)
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
    for n in eachindex(ions), m in eachindex(T.lasers)
        if m ≡ 1
            rn = length(ions) - n + 1
            ts = subleveltransitions(ions[n])
            C = sum([
                i * real.(sigma(ions[n], reverse(ts[i])...).data) for i in 1:length(ts)
            ])
            if length(ions) == 1
                K = C
            else
                K = kron(ion_arrays[1:(n - 1)]..., C, ion_arrays[(n + 1):length(ions)]...)
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
            for (i, mode) in enumerate(modes)
                push!(D_arrays, real.(displace(mode, ηlist(0)[i]).data))
            end
        end

        # iterate over ion-laser transitions
        for (ti, tr) in enumerate(ts)
            Δ, Ω = Δm[rn, m][ti], Ωm[rn, m][ti]
            Δ_2π = Δ / 2π
            typeof(Ω) <: Number && continue  # e.g. the laser doesn't shine on this ion
            locs = view(ion_idxs, ((ti - 1) * ion_reps + 1):(ti * ion_reps))
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