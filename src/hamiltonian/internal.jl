using SparseArrays: rowvals, nzrange, spzeros, spdiagm, findnz, rowvals, nonzeros
using FunctionWrappers: FunctionWrapper
using PolynomialRoots: roots
using QuantumOptics: SparseOperator, embed

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
    ions = T.configuration.ions
    vms = get_vibrational_modes(T.configuration)
    lasers = T.lasers
    (N, M, L) = map(x -> length(x), [ions, lasers, vms])
    ηnml = Array{Any}(undef, N, M, L)
    for n in 1:N, m in 1:M, l in 1:L
        δν = vms[l].δν
        ν = vms[l].ν
        eta = get_η(vms[l], lasers[m], ions[n], scaled = true)
        if eta == 0
            ηnml[n, m, L - l + 1] = 0
        else
            ηnml[n, m, L - l + 1] =
                FunctionWrapper{Float64, Tuple{Float64}}(t -> eta / √(ν + δν(t)))
        end
    end
    return ηnml
end

# Returns an array of vectors. The rows and columns of the array refer to ions and lasers,
# respectively. For each row/column we have a vector of detunings from the laser frequency
# for each ion transition. We need to separate this calculation from _Ωmatrix to implement
# RWA easily.
function _Δmatrix(T, timescale)
    ions = T.configuration.ions
    lasers = T.lasers
    (N, M) = length(ions), length(lasers)
    B = T.B
    ∇B = T.∇B
    Δnmkj = Array{Vector}(undef, N, M)
    for n in 1:N, m in 1:M
        Btot = B + ∇B * ionposition(ions[n])
        v = Vector{Float64}(undef, 0)
        for transition in subleveltransitions(ions[n])
            ωa = transitionfrequency(ions[n], transition, B = Btot)
            push!(v, 2π * timescale * ((c / lasers[m].λ) + lasers[m].Δ - ωa))
        end
        Δnmkj[n, m] = v
    end
    return Δnmkj
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
        transitions = subleveltransitions(ions[n])
        s_indx = findall(x -> x[1] == n, lasers[m].pointing)
        if length(s_indx) == 0
            Ωnmkj[n, m] = [0 for _ in 1:length(transitions)]
            continue
        else
            s = lasers[m].pointing[s_indx[1]][2]
        end
        v = []
        for t in transitions
            Ω0 =
                2π *
                timescale *
                s *
                matrix_element(ions[n], t, 1.0, lasers[m].k, lasers[m].ϵ, T.Bhat) / 2.0
            if Ω0 == 0
                push!(v, 0)
            else
                push!(
                    v,
                    FunctionWrapper{ComplexF64, Tuple{Float64}}(
                        t -> Ω0 * E(t) * exp(-im * phase(t))
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
    @simd for i in 0:(L - 1)
        if i == 0
            @inbounds rowcol .+= indxs[L - i]
        else
            @inbounds rowcol .+= (indxs[L - i] .- 1) .* prod(view(dims, 1:i))
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
    for i in (m + 1):n
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
