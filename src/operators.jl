using QuantumOptics: projector, tensor, SparseOperator, DenseOperator
using LinearAlgebra: diagm
import QuantumOptics: displace, thermalstate, coherentthermalstate, fockstate


export create, destroy, number, displace, coherentstate, coherentthermalstate, fockstate,
       thermalstate, sigma, ion_state, ion_projector


#############################################################################################
# vibrational_mode operators
#############################################################################################

"""
    create(v::vibrational_mode)
returns the creation operator for `v` such that: `create(v) * v[i] = √(i+1) * v[i+1]`.
"""
create(v::vibrational_mode) = SparseOperator(v, diagm(-1 => sqrt.(1:v.N)))

"""
    destroy(v::vibrational_mode)
Returns the destruction operator for `v` such that: `destroy(v) * v[i] = √i * v[i-1]`.
"""
destroy(v::vibrational_mode) = create(v)'

"""
    number(v::vibrational_mode)
Returns the number operator for `v` such that:  `number(v) * v[i] = i * v[i]`.
"""
number(v::vibrational_mode) = SparseOperator(v, diagm(0 => 0:v.N))

"""
    displace(v::vibrational_mode, α::Number)
Returns the displacement operator ``D(α)`` corresponding to `v`.
"""
function displace(v::vibrational_mode, α::Number)
    D = zeros(ComplexF64, v.N+1, v.N+1)
    @inbounds for n in 1:v.N+1, m in 1:v.N+1
        D[n, m] = _Dnm(α, n, m)
    end
    DenseOperator(v, D)
end

"""
    thermalstate(v::vibrational_mode, n̄::Real)
Returns a thermal density matrix with ``⟨a^†a⟩ ≈ n̄``. Note: approximate because we are 
dealing with a finite dimensional Hilbert space that must be normalized.
"""
function thermalstate(v::vibrational_mode, n̄::Real)
    if n̄ == 0
        return v[0] ⊗ v[0]'
    end
    T = 1 / (log((1 / n̄) + 1))
    H = create(v) * destroy(v)
    thermalstate(H, T)
end

"""
    coherentstate(v::vibrational_mode, α::Number)
Returns a coherent state on `v` with complex amplitude ``α``.
"""
function coherentstate(v::vibrational_mode, α::Number)
    # this implementation is the same as in QuantumOptics.jl, but there the function is 
    # restricted to v::FockBasis, so we must reimplement here
    k = zeros(ComplexF64, v.N+1)
    k[1] = exp(-abs2(α) / 2)
    @inbounds for n=1:v.N
        k[n+1] = k[n] * α / √n
    end
    k
end

"""
    coherentthermalstate(v::vibrational_mode, n̄::Real, α::Number)
Returns a displaced thermal state for `v`. The mean occupation of the thermal state is `n̄` 
, and `α` is the complex amplitude of the displacement.
"""
function coherentthermalstate(v::vibrational_mode, n̄::Real, α::Number)
    d = displace(v, α)
    d * thermalstate(v, n̄) * d'
end

"""
    fockstate(v::vibrational_mode, N::Int)
Returns the fockstate ``|N⟩`` on `v`. 
"""
fockstate(v::vibrational_mode, N::Int) = v[N]

#############################################################################################
# ion operators
#############################################################################################

"""
    sigma(ion::Ion, ψ1::String, ψ2::String)
<br>Creates ``|ψ1\\rangle\\langle ψ2|``, where ``|ψ_i\\rangle`` corresponds to the state
returned by `ion[ψᵢ]`.
"""
sigma(ion::Ion, ψ1::String, ψ2::String) = projector(ion[ψ1], dagger(ion[ψ2]))

"""
    ion_state(T::Trap, states...)
<br>i.e. `ion_state(T, "S-1/2", "S-1/2")` returns `ion1["S-1/2"] ⊗ ion2["S-1/2"]`
"""
function ion_state(T::trap, states...)
    N = length(T.configuration.ions)
    @assert N == length(states) "wrong number of states"
    tensor([T.configuration.ions[i][states[i]] for i in 1:N])
end

"""
    ion_projector(T::trap, states...) 
<br>i.e. `ion_projector(T, "S-1/2", "S-1/2")` returns:
``|S-1/2\\rangle\\langle S-1/2| ⊗ |S-1/2\\rangle\\langle S-1/2| ⊗ 𝐼``, 
where the identity ``𝐼`` is over the tensor product of all vibrational states.
"""
function ion_projector(T::trap, states...)
    N = length(T.configuration.ions)
    @assert N == length(states) "wrong number of states"
    modes = collect(Iterators.flatten(T.configuration.vibrational_modes))
    state = projector(T.configuration.ions[1][states[1]])
    for i in 2:N
        state = state ⊗ projector(T.configuration.ions[i][states[i]])
    end
    for mode in modes
        state = state ⊗ one(mode)
    end
    state    
end