using QuantumOptics
using Test, IonSim
using IonSim.analytical
using Suppressor
using Unitful

@suppress_err begin
    @testset "Dynamics -- Rabi" begin
        # carrier
        C = Ca40([("S1/2", -1 / 2), ("D5/2", -1 / 2)])
        L = Laser()
        chain = LinearChain(
            ions = [C],
            com_frequencies = (x = 3e6u"1/s", y = 3e6u"1/s", z = 1e6u"1/s"),
            vibrational_modes = (; z = [1])
        )
        T = Trap(configuration = chain, B = 4e-4u"T", Bhat = ẑ, δB = 0u"T", lasers = [L])
        L.λ = transitionwavelength(C, (("S1/2", -1 / 2), ("D5/2", -1 / 2)), T)
        mode = T.configuration.vibrational_modes.z[1]
        L.k = (x̂ + ẑ) / √2
        L.ϵ = (x̂ - ẑ) / √2 # I have no idea how this worked before; previously this was set to ŷ which gives zero coupling strength to Δm=0
        Ω = (rand() + 0.1) * 1e6
        Efield_from_rabi_frequency!(Ω*1u"1/s", T, 1, 1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)))

        h = hamiltonian(T)
        tspan = 0:1e-3:4
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[0],
            h
        )
        ex = expect(ionprojector(T, ("D5/2", -1 / 2)), sol)

        @test isapprox(real.(ex), @.(sin(2π * Ω / 2 * tout)^2), rtol = 1e-6)

        # add detuning
        Δ = (rand() + 0.5) * 10e5
        L.Δ = Δ*1u"1/s"
        h = hamiltonian(T)
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[0],
            h
        )
        ex = expect(ionprojector(T, ("D5/2", -1 / 2)), sol)
        @test isapprox(
            real.(ex),
            @.((Ω^2 / (Ω^2 + Δ^2)) * sin(2π * √(Ω^2 + Δ^2) / 2 * tout)^2),
            rtol = 1e-6
        )

        # add detuning using ion's stark_shift field
        L.Δ = 0u"1/s"
        C.stark_shift[("S1/2", -1 / 2)] = -Δ*1u"1/s" / 2
        C.stark_shift[("D5/2", -1 / 2)] = Δ*1u"1/s" / 2
        h = hamiltonian(T)
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[0],
            h
        )
        ex = expect(ionprojector(T, ("D5/2", -1 / 2)), sol)
        @test isapprox(
            real.(ex),
            @.((Ω^2 / (Ω^2 + Δ^2)) * sin(2π * √(Ω^2 + Δ^2) / 2 * tout)^2),
            rtol = 1e-6
        )

        # hot carrier
        zero_stark_shift!(C)
        mode.N = 100
        ψi_ion = C[("S1/2", -1 / 2)] ⊗ C[("S1/2", -1 / 2)]'
        n̄ = rand(1:20)
        ψi_mode = thermalstate(mode, n̄)
        ψi = ψi_ion ⊗ ψi_mode
        tspan = 0:1e-3:4
        h = hamiltonian(T, lamb_dicke_order = 0)
        @time tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h)
        ex = expect(ionprojector(T, ("D5/2", -1 / 2)), sol)
        η = get_η(mode, L, C)
        @test isapprox(real.(ex), rabi_flop(tout, Ω / 2, η, n̄), rtol = 1e-6)

        # sideband transitions

        ## under an RWA, RSB transitions, when starting in motional ground state,
        ## should be suppressed
        L.Δ = -mode.ν
        h = hamiltonian(T, rwa_cutoff = 1e5)
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[0],
            h
        )
        ex = expect(ionprojector(T, ("D5/2", -1 / 2)), sol)
        @test sum(ex) == 0

        ## a BSB should have a coupling strength ηΩ√(n+1)
        ## a RSB should have a coupling strength ηΩ√n
        L.Δ = -mode.ν
        h = hamiltonian(T, rwa_cutoff = 1e5)
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[1],
            h
        )
        ex = expect(ionprojector(T, ("D5/2", -1 / 2)), sol)
        @test isapprox(real.(ex), @.(sin(2π * η * Ω / 2 * tout)^2), rtol = 1e-6)

        L.Δ = mode.ν
        h = hamiltonian(T, rwa_cutoff = 1e5)
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[1],
            h
        )
        ex = expect(ionprojector(T, ("D5/2", -1 / 2)), sol)
        @test isapprox(real.(ex), @.(sin(2π * sqrt(2) * η * Ω / 2 * tout)^2), rtol = 1e-6)
    end

    @testset "Dynamics -- time-dependent-parameters" begin
        # could use some better tests

        # use ϕ(t) to set laser frequency to resonance
        C = Ca40([("S1/2", -1 / 2), ("D5/2", -1 / 2)])
        L = Laser()
        chain = LinearChain(
            ions = [C],
            com_frequencies = (x = 3e6u"1/s", y = 3e6u"1/s", z = 1e6u"1/s"),
            vibrational_modes = (; z = [1])
        )
        T = Trap(configuration = chain, B = 4e-4u"T", Bhat = ẑ, δB = 0u"T", lasers = [L])
        L.λ = transitionwavelength(C, (("S1/2", -1 / 2), ("D5/2", -1 / 2)), T)
        mode = T.configuration.vibrational_modes.z[1]
        L.k = (x̂ + ẑ) / √2
        L.ϵ = (x̂ - ẑ) / √2 # I have no idea how this worked before; previously this was set to ŷ which gives zero coupling strength to Δm=0
        Δf =
            transitionfrequency(1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)), T) -
            transitionfrequency(1, ("S1/2", "D5/2"), T)
        L.ϕ = t -> Δf * t
        Ω = (rand() + 0.1) * 1e6
        Efield_from_rabi_frequency!(Ω*1u"1/s", T, 1, 1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)))

        h = hamiltonian(T)
        tspan = 0:1e-3:4
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[0],
            h
        )
        ex = expect(ionprojector(T, ("D5/2", -1 / 2)), sol)

        @test isapprox(real.(ex), @.(sin(2π * Ω / 2 * tout)^2), rtol = 1e-6)

        # set Ω(t) to a step function
        L.Δ = 0u"1/s"
        L.ϕ = 0
        E = Efield_from_pi_time(Ω, T, 1, 1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)))
        L.E = t -> t < 1 ? 0u"V/m" : E
        h = hamiltonian(T)
        tspan = 0:1e-3:3
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[0],
            h
        )
        ex = expect(ionprojector(T, ("D5/2", -1 / 2)), sol)

        delayed_sin(t) = t < 1 ? 0 : sin(2π * Ω / 2 * (t - 1))
        @test isapprox(real.(ex), @.(delayed_sin(tout)^2), rtol = 1e-6)

        # check that δν(t) shifts sideband frequencies appropriately
        tspan = 0:1e-3:30
        mode.δν = t -> 20e3 * 1e-6
        L.Δ = mode.ν + 20e3
        h = hamiltonian(T, rwa_cutoff = 1e5)
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[1],
            h
        )
        ex = expect(ionprojector(T, ("D5/2", -1 / 2)), sol)
        η = get_η(mode, L, C)
        # η -> η(t) also
        @test isapprox(
            real.(ex),
            @.(sin(2π * sqrt(2) * (η / sqrt(1.02)) * Ω / 2 * tout)^2),
            rtol = 1e-6
        )

        # δB(t)
        ## add test
    end

    @testset "molmer-sorensen" begin
        # slow gate
        C = Ca40([("S1/2", -1 / 2), ("D5/2", -1 / 2)])
        L1 = Laser()
        L2 = Laser()
        chain = LinearChain(
            ions = [C, C],
            com_frequencies = (x = 3e6u"1/s", y = 3e6u"1/s", z = 1e6u"1/s"),
            vibrational_modes = (; z = [1])
        )
        T = Trap(configuration = chain, B = 4e-4u"T", Bhat = (x̂ + ẑ) / √2, lasers = [L1, L2])
        L1.λ = transitionwavelength(C, (("S1/2", -1 / 2), ("D5/2", -1 / 2)), T)
        L2.λ = transitionwavelength(C, (("S1/2", -1 / 2), ("D5/2", -1 / 2)), T)
        mode = T.configuration.vibrational_modes.z[1]
        # Set the laser parameters
        ϵ = 40e3
        d = 80  # corrects for AC stark shift from single-photon coupling to sidebands
        L1.Δ = mode.ν + ϵ - d
        L1.k = ẑ
        L1.ϵ = x̂

        L2.Δ = mode.ν - ϵ + d
        L2.k = ẑ
        L2.ϵ = x̂

        mode.N = 15
        η = abs(get_η(mode, L1, C))
        Ω = √(1e3 * ϵ) / η  # This will give a 1kHz MS strength, since coupling goes like (ηΩ)^2/ϵ

        Efield_from_rabi_frequency!(Ω*1u"1/s", T, 1, 1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)))
        Efield_from_rabi_frequency!(Ω*1u"1/s", T, 2, 1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)))
        ψi = ionstate(T, ("S1/2", -1 / 2), ("S1/2", -1 / 2)) ⊗ mode[0]  # initial state
        h = hamiltonian(T, rwa_cutoff = 5e5)
        tspan = 0:0.25:1000
        @time tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h)
        SS = expect(ionprojector(T, ("S1/2", -1 / 2), ("S1/2", -1 / 2)), sol)
        DD = expect(ionprojector(T, ("D5/2", -1 / 2), ("D5/2", -1 / 2)), sol)
        ex = analytical.two_ion_ms(tspan, 1e-6Ω, 1e-6mode.ν, 1e-6mode.ν + 1e-6ϵ, η, 0)
        @test isapprox(ex[1], SS, rtol = 1e-2)
        @test isapprox(ex[2], DD, rtol = 1e-2)

        # add checks for fast and hot gates
    end
end  # end suppress
