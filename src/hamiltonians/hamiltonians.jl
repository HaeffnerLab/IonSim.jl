export Hamiltonian

mutable struct Hamiltonian
    chamber::Chamber  # Pull user-defined params from here

    bare_hamiltonian
    interacting_hamiltonian
    rotatingframe
    rotatingframe_hamiltonian


    timescale::Real  # Scaling factor to be applied to all instances of time
    lambdickeorder::Union{Vector{Int}, Int}  # Only terms order η^lambdickeorder are retained
    rwacutoff::Real  # Time-dep terms with frequencies > rwacutoff are discarded
    displacement::String  # Approximation technique to be applied to displacement operators
    timedependenteta::Bool # Should η's reflect time-dependence of trap frequencies?
end

function Hamiltonian(chamber::Chamber, 
        bare_hamiltonian, interacting_hamiltonian, 
        rotating_frame, rotating_frame_hamiltonian
    )
    Hamiltonian(
        chamber,
        bare_hamiltonian, interacting_hamiltonian, 
        rotating_frame, rotating_frame_hamiltonian,
        1, 1, Inf, "truncated", false
    )
end

