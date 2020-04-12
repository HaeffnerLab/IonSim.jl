using QuantumOptics: projector, tensor, dagger, create, destroy
import QuantumOptics: thermalstate, coherentstate, coherentthermalstate


export sigma, ion_state, ion_projector


"""
    thermalstate(v::vibrational_mode, n̄::Real)
<br>returns a thermal density matrix with ``\\langle aa^{\\dagger}\\rangle = `` n
"""
thermalstate(v::vibrational_mode, n̄::Real) = thermalstate(v.basis, n̄)
    
function thermalstate(b::Basis, n̄::Real)
    if n̄ == 0
        return fockstate(b, 0) ⊗ dagger(fockstate(b, 0))
    end
    T = 1 / (log((1 / n̄) + 1))
    H = create(b) * destroy(b)
    thermalstate(H, T)
end

coherentstate(v::vibrational_mode, α::Number) = coherentstate(v.basis, α)

function coherentthermalstate(b::Basis, α)
    return
end

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
        state = state ⊗ one(mode.basis)
    end
    state    
end