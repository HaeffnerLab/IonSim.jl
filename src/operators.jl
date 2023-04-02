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
    sigma,
    ionprojector,
    ionstate

#############################################################################################
# VibrationalMode operators
#############################################################################################

"""
    create(v::VibrationalMode)
returns the creation operator for `v` such that: `create(v) * v[i] = √(i+1) * v[i+1]`.
"""
create(v::VibrationalMode) = SparseOperator(v, diagm(-1 => sqrt.(1:(modecutoff(v)))))

"""
    destroy(v::VibrationalMode)
Returns the destruction operator for `v` such that: `destroy(v) * v[i] = √i * v[i-1]`.
"""
destroy(v::VibrationalMode) = create(v)'

"""
    number(v::VibrationalMode)
Returns the number operator for `v` such that:  `number(v) * v[i] = i * v[i]`.
"""
number(v::VibrationalMode) = SparseOperator(v, diagm(0 => 0:(modecutoff(v))))

"""
    displace(v::VibrationalMode, α::Number; method="truncated")
Returns the displacement operator ``D(α)`` corresponding to `v`.

If `method="truncated"` (default), the matrix elements are computed according to
``D(α) = exp[αa^† - α^*a]`` where ``a`` and ``a^†`` live in a truncated Hilbert space of
dimension `modecutoff(v)+1`.
Otherwise if `method="analytic"`, the matrix elements are computed assuming an
infinite-dimension Hilbert space. In general, this option will not return a unitary operator.
"""
function displace(v::VibrationalMode, α::Number; method="truncated")
    # @assert v.N ≥ abs(α) "`α` must be less than `v.N`"
    # Above line commented out to allow for Hamiltonian construction even if vibrational mode N = 0.
    # May want to think of a different way to perform this check in the future.
    @assert method in ["truncated", "analytic"] "method ∉ [truncated, analytic]"
    D = zeros(ComplexF64, modecutoff(v) + 1, modecutoff(v) + 1)
    if α == 0
        return one(v)
    elseif method ≡ "analytic"
        @inbounds begin
            @simd for n in 1:(modecutoff(v)+1)
                @simd for m in 1:(modecutoff(v)+1)
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
function thermalstate(v::VibrationalMode, n̄::Real; method="truncated")
    @assert modecutoff(v) ≥ n̄ "`n̄` must be less than `modecutoff(v)`"
    @assert method in ["truncated", "analytic"] "method ∉ [truncated, analytic]"
    if n̄ == 0
        return v[0] ⊗ v[0]'
    elseif method ≡ "truncated"
        d = [(n̄ / (n̄ + 1))^i for i in 0:(modecutoff(v))]
        return DenseOperator(v, diagm(0 => d) ./ sum(d))
    elseif method ≡ "analytic"
        return DenseOperator(
            v,
            diagm(0 => [(n̄ / (n̄ + 1))^i / (n̄ + 1) for i in 0:(modecutoff(v))])
        )
    end
end

"""
    coherentstate(v::VibrationalMode, α::Number)
Returns a coherent state on `v` with complex amplitude ``α``.
"""
function coherentstate(v::VibrationalMode, α::Number)
    # this implementation is the same as in QuantumOptics.jl, but there the function is
    # restricted to v::FockBasis, so we must reimplement here
    @assert modecutoff(v) ≥ abs(α) "`α` must be less than `modecutoff(v)`"
    k = zeros(ComplexF64, modecutoff(v) + 1)
    k[1] = exp(-abs2(α) / 2)
    @inbounds for n in 1:(modecutoff(v))
        k[n+1] = k[n] * α / √n
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
function coherentthermalstate(v::VibrationalMode, n̄::Real, α::Number; method="truncated")
    @assert (modecutoff(v) ≥ n̄ && modecutoff(v) ≥ abs(α)) "`n̄`, `α` must be less than `modecutoff(v)`"
    @assert method in ["truncated", "analytic"] "method ∉ [truncated, analytic]"
    if method ≡ "truncated"
        d = displace(v, α)
    elseif method ≡ "analytic"
        d = displace(v, α, method="analytic")
    end
    return d * thermalstate(v, n̄) * d'
end

"""
    fockstate(v::VibrationalMode, N::Int)
Returns the fockstate ``|N⟩`` on `v`.
"""
fockstate(v::VibrationalMode, N::Int) = v[N]

#############################################################################################
# Ion operators
#############################################################################################

"""
    ionstate(ion::Ion, sublevel)
Retuns the ket corresponding to the `Ion` being in state ``|sublevel⟩``. Options:
sublevel <: Tuple{String,Real}: Specifies full sublevel name
sublevel <: String: Specifies sublevel alias
sublevel <: Int: Returns the `sublevel`th eigenstate
"""
function ionstate(ion::Ion, sublevel::Tuple{String, Real})
    validatesublevel(ion, sublevel)
    i = findall(sublevels(ion) .== [sublevel])[1]
    return basisstate(ion, i)
end
ionstate(ion::Ion, sublevel::Int) = basisstate(ion, sublevel)

"""
    ionstate(object::Union{IonTrap, Chamber}, sublevels)
If `N = length(ions(object))`, returns N-dimensional ket corresponding to the ions being in
the state ``|sublevel₁⟩⊗|sublevel₂⟩⊗...⊗|sublevel\\_N⟩``.

`sublevels` must be an length-`N` Vector, with each element specifying its corresponding
ion's sublevel, using the same syntax as in `ionstate(ion::Ion, sublevel)`.
"""
function ionstate(iontrap::IonTrap, states::Vector)
    allions = ions(iontrap)
    L = length(allions)
    @assert L ≡ length(states) "wrong number of states"
    return tensor([ionstate(allions[i], states[i]) for i in 1:L])
end
ionstate(chamber::Chamber, states::Vector) = ionstate(iontrap(chamber), states)

"""
    sigma(ion::Ion, ψ1::sublevel[, ψ2::sublevel])
Returns ``|ψ₁⟩⟨ψ₂|``, where ``|ψᵢ⟩`` corresponds to the state returned by `ion[ψᵢ]`.

If ψ2 is not given, then the projector ``|ψ₁⟩⟨ψ₁|`` is returned.
"""
sigma(ion::Ion, ψ1::T, ψ2::T) where {T <: Union{Tuple{String, Real}, String, Int}} =
    sparse(projector(ion[ψ1], dagger(ion[ψ2])))
sigma(ion::Ion, ψ1::Union{Tuple{String, Real}, String, Int}) = sigma(ion, ψ1, ψ1)

"""
    ionprojector(obj, sublevels...; only_ions=false)

If `obj<:IonTrap` this will return ``|ψ₁⟩⟨ψ₁|⊗...⊗|ψ\\_N⟩⟨ψ\\_N|⊗𝟙``
where ``|ψᵢ⟩`` = `obj.ions[i][sublevels[i]]` and the identity operator ``𝟙`` is over all of 
the motional modes considered in `obj`.

If `only_ions=true`, then the projector is defined only over the ion subspace.

If instead `obj<:Chamber`, then this is the same as `obj = Chamber.iontrap`.
"""
function ionprojector(
    IC::IonTrap,
    sublevels::Union{Tuple{String, Real}, String, Int}...;
    only_ions=false
)
    allions = ions(IC)
    L = length(allions)
    @assert L ≡ length(sublevels) "wrong number of sublevels"
    allmodes = modes(IC)
    observable = tensor([projector(allions[i][sublevels[i]]) for i in 1:L]...)
    if !only_ions
        for mode in allmodes
            observable = observable ⊗ one(mode)
        end
    end
    return observable
end

function ionprojector(
    T::Chamber,
    sublevels::Union{Tuple{String, Real}, String, Int}...;
    only_ions=false
)
    return ionprojector(iontrap(T), sublevels..., only_ions=only_ions)
end

#############################################################################################
# internal functions
#############################################################################################

# computes iⁿ(-i)ᵐ * (s! / ((s+1) * √(m!n!)))
function _pf(s::Int, n::Int, m::Int)
    n -= 1
    m -= 1
    s -= 1
    @assert n <= s && m <= s
    val = 1.0 / (s + 1)
    for i in 0:(s-2)
        if (m - i > 0) && (n - i > 0)
            val *= (s - i) / (√((m - i) * (n - i)))
        elseif m - i > 0
            val *= (s - i) / (√(m - i))
        elseif n - i > 0
            val *= (s - i) / (√(n - i))
        else
            val *= (s - i)
        end
    end
    return (-1im)^n * 1im^m * val
end

# computes the coefficients for the 'probabilist's' Hermite polynomial of order n
function _He(n::Int)
    a = zeros(Float64, n + 2, n + 2)
    a[1, 1] = 1
    a[2, 1] = 0
    a[2, 2] = 1
    for i in 2:(n+1), j in 1:(n+1)
        if j ≡ 1
            a[i+1, j] = -(i - 1) * a[i-1, j]
        else
            a[i+1, j] = a[i, j-1] - (i - 1) * a[i-1, j]
        end
    end
    return [a[n+1, k+1] for k in 0:n]
end

# computes He_n(x) (nth order Hermite polynomial)
function _fHe(x::Real, n::Int)
    n -= 1
    He = 1.0, x
    if n < 2
        return He[n+1]
    end
    for i in 2:n
        He = He[2], x * He[2] - (i - 1) * He[1]
    end
    return He[2]
end

# computes the matrix elements ⟨m|Dˢ(α)|n⟩ for the truncated displacement operator Dˢ(α)
# which exists in a Hilbert space of dimension s
function _Dtrunc(Ω, Δ, η, ν, rs, s, n, prefactor, timescale, L, t)
    d = complex(1, 0)
    for i in 1:L
        val = 0.0
        Δn = n[1][i] - n[2][i]
        for r in rs[i]
            val +=
                exp(im * r * abs(η[i])) * _fHe(r, n[2][i]) * _fHe(r, n[1][i]) /
                _fHe(r, s[i])^2
        end
        d *= (
            exp(im * Δn * (2π * ν[i] * timescale * t + π / 2 + π * (sign(η[i] < 0)))) *
            val *
            prefactor[i]
        )
    end
    g = Ω * exp(-1im * t * Δ)
    return g * d, g * conj(d)
end

# associated Laguerre polynomial
function _alaguerre(x::Real, n::Int, k::Int)
    L = 1.0, -x + k + 1
    if n < 2
        return L[n+1]
    end
    for i in 2:n
        L = L[2], ((k + 2i - 1 - x) * L[2] - (k + i - 1) * L[1]) / i
    end
    return L[2]
end

# matrix elements of the displacement operator in the Fock Basis, assuming an
# infinite-dimensional Hilbert space. https://doi.org/10.1103/PhysRev.177.1857
function _Dnm(ξ::Number, n::Int, m::Int)
    if n < m
        return (-1)^isodd(abs(n - m)) * conj(_Dnm(ξ, m, n))
    end
    n -= 1
    m -= 1
    s = 1.0
    for i in (m+1):n
        s *= i
    end
    ret = sqrt(1 / s) * ξ^(n - m) * exp(-abs2(ξ) / 2.0) * _alaguerre(abs2(ξ), m, n - m)
    if isnan(ret)
        return 1.0 * (n == m)
    end
    return ret
end
