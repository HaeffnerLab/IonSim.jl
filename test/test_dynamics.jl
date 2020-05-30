using QuantumOptics
using Test, IonSim
using IonSim.analytical


@testset "rabi" begin
    # carrier
    C = ca40(selected_level_structure=["S-1/2", "D-1/2"])
    L = laser()
    chain = linearchain(
            ions=[C], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[], y=[], z=[1])
            )
    T = trap(configuration=chain, B=4e-4, Bhat=ẑ, δB=0, lasers=[L])
    mode = T.configuration.vibrational_modes.z[1]
    L.k = (x̂ + ẑ)/√2 
    L.ϵ = ŷ
    Δf = transition_frequency(T, 1, ("S-1/2", "D-1/2"))
    L.Δ = Δf
    Efield_from_pi_time!(2e-6, T, 1, 1, ("S-1/2", "D-1/2"))
    
    h = hamiltonian(T)
    tspan = 0:0.01:4
    tout, sol = timeevolution.schroedinger_dynamic(
            tspan, ion_state(T, "S-1/2") ⊗ fockstate(mode.basis, 0), h
        )
    ex = expect(ion_projector(T, "D-1/2"), sol)
    
    @test isapprox(real.(ex), @.(sin(2π * 0.25/2 * tout)^2), rtol=1e-2)

    # add detuning
    L.Δ = Δf + 2.5e5
    h = hamiltonian(T)
    tout, sol = timeevolution.schroedinger_dynamic(
            tspan, ion_state(T, "S-1/2") ⊗ fockstate(mode.basis, 0), h
        )
    ex = expect(ion_projector(T, "D-1/2"), sol)
    @test isapprox(real.(ex), @.((1/2) * sin(2π * √2 * 0.25/2 * tout)^2), rtol=1e-2)

    # add detuning using ion's stark_shift field
    L.Δ = Δf
    C.stark_shift["S-1/2"] = -1.25e5
    C.stark_shift["D-1/2"] = 1.25e5
    h = hamiltonian(T)
    tout, sol = timeevolution.schroedinger_dynamic(
            tspan, ion_state(T, "S-1/2") ⊗ fockstate(mode.basis, 0), h
        )
    ex = expect(ion_projector(T, "D-1/2"), sol)
    @test isapprox(real.(ex), @.((1/2) * sin(2π * √2 * 0.25/2 * tout)^2), rtol=1e-2)

    # hot carrier
    zero_stark_shift(C)
    mode.N = 50
    ψi_ion = C["S-1/2"] ⊗ C["S-1/2"]'
    ψi_mode = thermalstate(mode, 10)
    ψi = ψi_ion ⊗ ψi_mode
    tspan = 0:0.1:4
    h = hamiltonian(T, lamb_dicke_order=0)
    @time tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h)
    ex = expect(ion_projector(T, "D-1/2"), sol)
    η = get_η(mode, L, C)
    @test isapprox(real.(ex), rabi_flop(tout, 0.25/2, η, 10), rtol=1e-2)

    # sideband transitions
    
    ## under an RWA, RSB transitions, when starting in motional ground state, 
    ## should be suppressed
    L.Δ = Δf - mode.ν
    h = hamiltonian(T, rwa_cutoff=1e5)
    tout, sol = timeevolution.schroedinger_dynamic(
            tspan, ion_state(T, "S-1/2") ⊗ fockstate(mode.basis, 0), h
        )
    ex = expect(ion_projector(T, "D-1/2"), sol)
    @test sum(ex) == 0

    ## a BSB should have a coupling strength ηΩ√(n+1)
    ## a RSB should have a coupling strength ηΩ√n
    L.Δ = Δf - mode.ν
    h = hamiltonian(T, rwa_cutoff=1e5)
    tout, sol = timeevolution.schroedinger_dynamic(
            tspan, ion_state(T, "S-1/2") ⊗ fockstate(mode.basis, 1), h
        )
    ex = expect(ion_projector(T, "D-1/2"), sol)
    @test isapprox(real.(ex), @.(sin(2π * η * 0.25/2 * tout)^2), rtol=1e-2)

    L.Δ = Δf + mode.ν
    h = hamiltonian(T, rwa_cutoff=1e5)
    tout, sol = timeevolution.schroedinger_dynamic(
            tspan, ion_state(T, "S-1/2") ⊗ fockstate(mode.basis, 1), h
        )
    ex = expect(ion_projector(T, "D-1/2"), sol)
    @test isapprox(real.(ex), @.(sin(2π * sqrt(2) * η * 0.25/2 * tout)^2), rtol=1e-2)
end

@testset "time-dependent-parameters" begin
    # could use some better tests

    # use ϕ(t) to set laser frequency to resonance
    C = ca40(selected_level_structure=["S-1/2", "D-1/2"])
    L = laser()
    chain = linearchain(
            ions=[C], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[], y=[], z=[1])
            )
    T = trap(configuration=chain, B=4e-4, Bhat=ẑ, δB=0, lasers=[L])
    mode = T.configuration.vibrational_modes.z[1]
    L.k = (x̂ + ẑ)/√2 
    L.ϵ = ŷ
    Δf = transition_frequency(T, 1, ("S-1/2", "D-1/2"))
    L.ϕ = t -> Δf * t
    Efield_from_pi_time!(2e-6, T, 1, 1, ("S-1/2", "D-1/2"))
    
    h = hamiltonian(T)
    tspan = 0:0.01:4
    tout, sol = timeevolution.schroedinger_dynamic(
            tspan, ion_state(T, "S-1/2") ⊗ fockstate(mode.basis, 0), h
        )
    ex = expect(ion_projector(T, "D-1/2"), sol)
    
    @test isapprox(real.(ex), @.(sin(2π * 0.25/2 * tout)^2), rtol=1e-2)

    # set Ω(t) to a step function
    L.Δ = Δf
    L.ϕ = 0
    E = Efield_from_pi_time(2e-6, T, 1, 1, ("S-1/2", "D-1/2"))
    L.E = t -> t < 1 ? 0 : E
    h = hamiltonian(T)
    tspan = 0:0.01:3
    tout, sol = timeevolution.schroedinger_dynamic(
            tspan, ion_state(T, "S-1/2") ⊗ fockstate(mode.basis, 0), h
        )
    ex = expect(ion_projector(T, "D-1/2"), sol)
    
    delayed_sin(t) = t < 1 ? 0 : sin(2π * 0.25/2 * (t-1))
    @test isapprox(real.(ex), @.(delayed_sin(tout)^2), rtol=1e-2)

    # check that δν(t) shifts sideband frequencies appropriately
    tspan = 0:0.01:30
    mode.δν = t -> 20e3 * 1e-6
    L.Δ = Δf + mode.ν + 20e3
    h = hamiltonian(T, rwa_cutoff=1e5)
    tout, sol = timeevolution.schroedinger_dynamic(
            tspan, ion_state(T, "S-1/2") ⊗ fockstate(mode.basis, 1), h
        )
    ex = expect(ion_projector(T, "D-1/2"), sol)
    η = get_η(mode, L, C)
                                                    # η -> η(t) also
    @test isapprox(real.(ex), @.(sin(2π * sqrt(2) * (η / sqrt(1.02)) * 0.25/2 * tout)^2), rtol=1e-1)

    # δB(t)
    ## add test
end

@testset "molmer-sorensen" begin
    # slow gate
    C = ca40(selected_level_structure=["S-1/2", "D-1/2"])
    L1 = laser()
    L2 = laser() 
    chain = linearchain(ions=[C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[], y=[], z=[1]))
    T = trap(configuration=chain, B=4e-4, Bhat=(x̂ + ẑ)/√2, lasers=[L1, L2])
    mode = T.configuration.vibrational_modes.z[1]
    # Set the laser parameters
    ϵ = 40e3
    d = 80  # corrects for AC stark shift from single-photon coupling to sidebands
    Δf = transition_frequency(T, 1, ("S-1/2", "D-1/2"))
    L1.Δ = Δf + mode.ν + ϵ - d
    L1.k = ẑ
    L1.ϵ = x̂

    L2.Δ = Δf - mode.ν - ϵ + d
    L2.k = ẑ
    L2.ϵ = x̂

    mode.N = 5
    η = abs(get_η(mode, L1, C))
    Ω = √(1e3 * ϵ) / η  # This will give a 1kHz MS strength, since coupling goes like (ηΩ)^2/ϵ

    Efield_from_rabi_frequency!(Ω, T, 1, 1, ("S-1/2", "D-1/2"))
    Efield_from_rabi_frequency!(Ω, T, 2, 1, ("S-1/2", "D-1/2"))
    ψi = ion_state(T, "S-1/2", "S-1/2") ⊗ fockstate(mode.basis, 0);  # initial state
    h = hamiltonian(T, rwa_cutoff=5e5)
    tspan = 0:0.25:1000
    @time tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h);
    SS = expect(ion_projector(T, "S-1/2", "S-1/2"), sol)
    DD = expect(ion_projector(T, "D-1/2", "D-1/2"), sol)
    ex = analytical.two_ion_ms(tspan, 1e-6Ω, 1e-6mode.ν, 1e-6mode.ν + 1e-6ϵ, η, 0)
    @test isapprox(ex[1], SS, rtol=1e-2)
    @test isapprox(ex[2], DD, rtol=1e-2)

    # add checks for fast and hot gates
end