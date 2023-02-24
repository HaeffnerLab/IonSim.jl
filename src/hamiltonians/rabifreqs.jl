export Ω0
export η
"""
This function computes the time-dependent Rabi frequency for a specific
transition in a specific ion driven by a specific laser. Note that the time dependence here is due to the fact that the laser's intensity may be a function of time. The complex exponential will be handled by a separate function.

The Rabi frequency Ω0 is calculated according to the equation

        Ω0 = 2πs(μ⃗⋅E⃗/ħ)

Here, μ⃗⋅E⃗ is the matrix element

        eE⟨j|r̂|i⟩

where e is the elementary charge, E is the magnitude of the laser's electric field, the transition in question is |i⟩→|j⟩ and r̂ is the position operator. The quantity μ, which is defined according to

        μ = e⟨j|r̂|i⟩

is evaluated by the matrixelement function. We also have that the electric field is related to the intensity I by E = √I, so that our expression reduces to

        Ω0(t) = 2πsμ√(I(t))

Here, s is a scale factor between 0 and 1 which may be set by the user to simulate an ion not seeing the full electric field of the laser.

In the lab frame, the interaction term has the form

        Ω0(t)[e^{iνt} + e^{-iνt}][|i⟩⟨j| + |j⟩⟨i|]

The function Ω0 returns Ω0 only, not the phasors.


"""

function Ω0(
        chamber,
        timescale,
        ion_idx,
        laser_idx,
        transition_idx
)
    ion = ions(chamber)[ion_idx]
    laser = lasers(chamber)[laser_idx]
    transition = subleveltransitions(ion)[transition_idx]
    I = intensity(laser)
    ϕ = phase(laser)
    s_idx = findall(x -> x[1] == ion_idx, pointing(laser))
    if length(s_idx) == 0
        return 0
    end
    s = pointing(laser)[s_idx[1]][2]

    Ω0 =
        2π *
        timescale *
        s *
        matrixelement(
                ion,
                transition,
                1.0,
                polarization(laser),
                wavevector(laser),
                bfield_unitvector(chamber)
        ) / 2.0
    if Ω0 == 0
        return 0
    end

    return FunctionWrapper{ComplexF64, Tuple{Float64}}(
                    t -> Ω0 * √I(t)
            )
end    

function η(
        chamber,
        ion_idx,
        laser_idx,
        mode_idx
    )
    
    ion = ions(chamber)[ion_idx]
    laser = lasers(chamber)[laser_idx]
    mode = modes(chamber)[mode_idx]
    
    δν = frequency_fluctuation(mode)
    ν = frequency(mode)
    eta = lambdicke(mode, ion, laser, scaled=true)
    if eta == 0
        return 0
    end
    return FunctionWrapper{Float64, Tuple{Float64}}(t -> eta / √(ν + δν(t)))
end
    
        