using QuantumOptics: projector, tensor, SparseOperator, DenseOperator, basisstate, Ket
using LinearAlgebra: diagm
import QuantumOptics: displace, thermalstate, coherentthermalstate, fockstate

export create,
    destroy,
    number,
    displace,
    coherentstate,
    coherentthermalstate,
    fockstate,
    thermalstate,
    get_η

"""
    get_η(V::VibrationalMode, L::Laser, I::Ion)
The Lamb-Dicke parameter:
``|k|cos(\\theta)\\sqrt{\\frac{\\hbar}{2m\\nu}}``
for a given vibrational mode, ion and laser.
"""
function get_η(V::VibrationalMode, L::Laser, I::Ion; scaled = false)
    @fastmath begin
        k = 2π / L.λ
        scaled ? ν = 1 : ν = V.ν
        x0 = √(ħ / (2 * mass(I) * 2π * ν))
        cosθ = ndot(L.k, V.axis)
        k * x0 * cosθ * V.mode_structure[ionnumber(I)]
    end
end

"""
    create(v::VibrationalMode)
returns the creation operator for `v` such that: `create(v) * v[i] = √(i+1) * v[i+1]`.
"""
create(v::VibrationalMode) = SparseOperator(v, diagm(-1 => sqrt.(1:(v.N))))

"""
    destroy(v::VibrationalMode)
Returns the destruction operator for `v` such that: `destroy(v) * v[i] = √i * v[i-1]`.
"""
destroy(v::VibrationalMode) = create(v)'

"""
    number(v::VibrationalMode)
Returns the number operator for `v` such that:  `number(v) * v[i] = i * v[i]`.
"""
number(v::VibrationalMode) = SparseOperator(v, diagm(0 => 0:(v.N)))

"""
    displace(v::VibrationalMode, α::Number; method="truncated")
Returns the displacement operator ``D(α)`` corresponding to `v`.

If `method="truncated"` (default), the matrix elements are computed according to
``D(α) = exp[αa^† - α^*a]`` where ``a`` and ``a^†`` live in a truncated Hilbert space of
dimension `v.N+1`.
Otherwise if `method="analytic"`, the matrix elements are computed assuming an
infinite-dimension Hilbert space. In general, this option will not return a unitary operator.
"""
function displace(v::VibrationalMode, α::Number; method = "truncated")
    # @assert v.N ≥ abs(α) "`α` must be less than `v.N`"
    # Above line commented out to allow for Hamiltonian construction even if vibrational mode N = 0.
    # May want to think of a different way to perform this check in the future.
    @assert method in ["truncated", "analytic"] "method ∉ [truncated, analytic]"
    D = zeros(ComplexF64, v.N + 1, v.N + 1)
    if α == 0
        return one(v)
    elseif method ≡ "analytic"
        @inbounds begin
            @simd for n in 1:(v.N + 1)
                @simd for m in 1:(v.N + 1)
                    D[n, m] = _Dnm(α, n, m)
                end
            end
        end
        return DenseOperator(v, D)
    elseif method ≡ "truncated"
        return exp(dense(α * create(v) - conj(α) * destroy(v)))
    end
end

"""
    thermalstate(v::VibrationalMode, n̄::Real; method="truncated")
Returns a thermal density matrix with ``⟨a^†a⟩ ≈ n̄``. Note: approximate because we are
dealing with a finite dimensional Hilbert space that must be normalized.

`method` can be set to either `"truncated"` (default) or `"analytic"`. In the former case,
the thermal density matrix is generated according to the formula:
``ρ_{th} = exp(-νa^†a/T) / Tr [exp(-νa^†a/T)]``. In the later case, the analytic formula,
assuming an infinite-dimensional Hilbert space, is used:
``[ρ_{th}]_{ij} = δ_{ij} \\frac{nⁱ}{(n+1)^{i+1}}.``
"""
function thermalstate(v::VibrationalMode, n̄::Real; method = "truncated")
    @assert v.N ≥ n̄ "`n̄` must be less than `v.N`"
    @assert method in ["truncated", "analytic"] "method ∉ [truncated, analytic]"
    if n̄ == 0
        return v[0] ⊗ v[0]'
    elseif method ≡ "truncated"
        d = [(n̄ / (n̄ + 1))^i for i in 0:(v.N)]
        return DenseOperator(v, diagm(0 => d) ./ sum(d))
    elseif method ≡ "analytic"
        return DenseOperator(v, diagm(0 => [(n̄ / (n̄ + 1))^i / (n̄ + 1) for i in 0:(v.N)]))
    end
end

"""
    coherentstate(v::VibrationalMode, α::Number)
Returns a coherent state on `v` with complex amplitude ``α``.
"""
function coherentstate(v::VibrationalMode, α::Number)
    # this implementation is the same as in QuantumOptics.jl, but there the function is
    # restricted to v::FockBasis, so we must reimplement here
    @assert v.N ≥ abs(α) "`α` must be less than `v.N`"
    k = zeros(ComplexF64, v.N + 1)
    k[1] = exp(-abs2(α) / 2)
    @inbounds for n in 1:(v.N)
        k[n + 1] = k[n] * α / √n
    end
    return Ket(v, k)
end

"""
    coherentthermalstate(v::VibrationalMode, n̄::Real, α::Number; method="truncated)
Returns a displaced thermal state for `v`, which is created by applying a displacement
operation to a thermal state. The mean occupation of the thermal state is `n̄` and `α` is the
complex amplitude of the displacement.

`method` can be either `"truncated"` or `"analytic"` and this argument determines how the
displacement operator is computed (see: [`displace`](@ref)) .
"""
function coherentthermalstate(v::VibrationalMode, n̄::Real, α::Number; method = "truncated")
    @assert (v.N ≥ n̄ && v.N ≥ abs(α)) "`n̄`, `α` must be less than `v.N`"
    @assert method in ["truncated", "analytic"] "method ∉ [truncated, analytic]"
    if method ≡ "truncated"
        d = displace(v, α)
    elseif method ≡ "analytic"
        d = displace(v, α, method = "analytic")
    end
    return d * thermalstate(v, n̄) * d'
end

"""
    fockstate(v::VibrationalMode, N::Int)
Returns the fockstate ``|N⟩`` on `v`.
"""
fockstate(v::VibrationalMode, N::Int) = v[N]
