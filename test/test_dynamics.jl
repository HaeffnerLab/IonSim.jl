using QuantumOptics
using Test, IonSim
using IonSim.analytical
using Suppressor

@suppress_err begin
    @testset "Dynamics -- Rabi" begin
        # carrier
        C = Ca40([("S1/2", -1 / 2), ("D5/2", -1 / 2)])
        L = Laser()
        chain = LinearChain(
            ions = [C],
            com_frequencies = (x = 3e6, y = 3e6, z = 1e6),
            selected_modes = (; z = [1])
        )
        T = Chamber(iontrap = chain, B = 4e-4, Bhat = ẑ, δB = 0, lasers = [L])
        L.λ = transitionwavelength(C, (("S1/2", -1 / 2), ("D5/2", -1 / 2)), T)
        mode = T.iontrap.selected_modes.z[1]
        L.k = (x̂ + ẑ) / √2
        L.ϵ = (x̂ - ẑ) / √2
        Ω = (rand() + 0.2) * 1e4 # Small so that sideband transitions are suppressed
        Ω00 = Ω * exp(-lambdicke(mode, L, C)^2 / 2) # Actual Rabi frequency of n=0 carrier Rabi oscillations
        Efield_from_rabi_frequency!(Ω, T, 1, 1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)))

        h = hamiltonian(T, timescale=1e-6)
        tspan = 0:1e-1:400
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[0],
            h
        )
        ex_ionsim_c0 = real.(expect(ionprojector(T, ("D5/2", -1 / 2)), sol))
        ex_analyt_c0 = @.(sin(2π * Ω00 / 2 * tout * 1e-6)^2)
        @test isapprox(ex_ionsim_c0, ex_analyt_c0, rtol = 1e-5)

        # This test serves to check for the presence of a bug -- The norm of the output seems to deviate from 1, moreso as the simulation goes on
        # If this test fails then the problem has likely been solved.
        # Find that the problem is worse if lamb_dicke_order=0
        @test maximum(abs.(1 .- norm.(sol))) > 1e-8

        # add detuning
        Δ = (rand() + 0.5) * 10e3
        L.Δ = Δ
        h = hamiltonian(T, timescale=1e-6)
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[0],
            h
        )
        ex_ionsim_d1 = real.(expect(ionprojector(T, ("D5/2", -1 / 2)), sol))
        ex_analyt_d1 =
            @.((Ω00^2 / (Ω00^2 + Δ^2)) * sin(2π * √(Ω00^2 + Δ^2) / 2 * tout * 1e-6)^2)
        @test isapprox(ex_ionsim_d1, ex_analyt_d1, rtol = 1e-4)

        # add detuning using ion's manual_shift
        L.Δ = 0
        C.manual_shift[("S1/2", -1 / 2)] = -Δ / 2
        C.manual_shift[("D5/2", -1 / 2)] = Δ / 2
        h = hamiltonian(T, timescale=1e-6)
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[0],
            h
        )
        ex_ionsim_d2 = real.(expect(ionprojector(T, ("D5/2", -1 / 2)), sol))
        ex_analyt_d2 =
            @.((Ω00^2 / (Ω00^2 + Δ^2)) * sin(2π * √(Ω00^2 + Δ^2) / 2 * tout * 1e-6)^2)
        @test isapprox(ex_ionsim_d2, ex_analyt_d2, rtol = 1e-4)

        # hot carrier
        # For this test, numerical result deviates from analytical quickly as nbar grows.
        # Keeping max at 10 ensures agreement to 10^-4; if max nbar is 20, find that need rtol no less than 10^-2
        zero_manual_shift!(C)
        mode.N = 100
        ψi_ion = C[("S1/2", -1 / 2)] ⊗ C[("S1/2", -1 / 2)]'
        n̄ = rand(1:10)
        ψi_mode = thermalstate(mode, n̄)
        ψi = ψi_ion ⊗ ψi_mode
        h = hamiltonian(T, timescale=1e-6, lamb_dicke_order = 0)
        tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h)
        η = lambdicke(mode, L, C)
        ex_ionsim_cn = real.(expect(ionprojector(T, ("D5/2", -1 / 2)), sol))
        ex_analyt_cn = analytical.rabi_flop(1e-6 * tout, Ω, η, n̄)
        @test isapprox(ex_ionsim_cn, ex_analyt_cn, rtol = 1e-2)

        # sideband transitions

        ## under an RWA, RSB transitions, when starting in motional ground state,
        ## should be suppressed
        L.Δ = -mode.ν
        h = hamiltonian(T, timescale=1e-6, rwa_cutoff = 1e3)
        tspan_sb = 0:1:2000
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan_sb,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[0],
            h
        )
        ex_ionsim_rsb0 = expect(ionprojector(T, ("D5/2", -1 / 2)), sol)
        @test sum(ex_ionsim_rsb0) == 0

        ## a BSB should have a coupling strength ηΩ√(n+1)
        ## a RSB should have a coupling strength ηΩ√n
        L.Δ = -mode.ν
        h = hamiltonian(T, timescale=1e-6, rwa_cutoff = 1e3)
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan_sb,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[1],
            h
        )
        ex_ionsim_rsb1 = real.(expect(ionprojector(T, ("D5/2", -1 / 2)), sol))
        ex_analyt_rsb1 = @.(sin(2π * η * Ω00 / 2 * tout * 1e-6)^2)
        @test isapprox(ex_ionsim_rsb1, ex_analyt_rsb1, rtol = 1e-5)

        L.Δ = mode.ν
        h = hamiltonian(T, timescale=1e-6, rwa_cutoff = 1e3)
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[1],
            h
        )
        ex_ionsim_bsb1 = real.(expect(ionprojector(T, ("D5/2", -1 / 2)), sol))
        ex_analyt_bsb1 = @.(sin(2π * sqrt(2) * η * Ω00 / 2 * tout * 1e-6)^2)
        @test isapprox(ex_ionsim_bsb1, ex_analyt_bsb1, rtol = 2e-2)
    end

    @testset "Dynamics -- time-dependent-parameters" begin
        # could use some better tests

        # use ϕ(t) to set laser frequency to resonance
        C = Ca40([("S1/2", -1 / 2), ("D5/2", -1 / 2)])
        L = Laser()
        chain = LinearChain(
            ions = [C],
            com_frequencies = (x = 3e6, y = 3e6, z = 1e6),
            selected_modes = (; z = [1])
        )
        T = Chamber(iontrap = chain, B = 4e-4, Bhat = ẑ, δB = 0, lasers = [L])
        L.λ = transitionwavelength(C, ("S1/2", "D5/2"), T)
        mode = T.iontrap.selected_modes.z[1]
        L.k = (x̂ + ẑ) / √2
        L.ϵ = (x̂ - ẑ) / √2 # I have no idea how this worked before; previously this was set to ŷ which gives zero coupling strength to Δm=0
        Δf =
            transitionfrequency(1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)), T) -
            transitionfrequency(1, ("S1/2", "D5/2"), T)
        L.ϕ = t -> 1e-6 * 2 * pi * Δf * t
        Ω = (rand() + 0.1) * 1e6
        Ω00 = Ω * exp(-lambdicke(mode, L, C)^2 / 2)
        Efield_from_rabi_frequency!(Ω, T, 1, 1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)))

        h = hamiltonian(T, timescale=1e-6, lamb_dicke_order = 0)
        tspan = 0:1e-3:4
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[0],
            h
        )
        ex_ionsim_tdphi = real.(expect(ionprojector(T, ("D5/2", -1 / 2)), sol))
        ex_analyt_tdphi = @.(sin(2π * 1e-6 * Ω00 / 2 * tout)^2)
        @test isapprox(ex_ionsim_tdphi, ex_analyt_tdphi, rtol = 1e-5)

        # set Ω(t) to a step function
        L.λ = transitionwavelength(1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)), T)
        L.ϕ = 0
        E = Efield_from_rabi_frequency(Ω, T, 1, 1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)))
        L.E = t -> t < 1 ? 0 : E
        h = hamiltonian(T, timescale=1e-6, lamb_dicke_order = 0)
        tspan = 0:1e-3:3
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[0],
            h
        )
        ex = expect(ionprojector(T, ("D5/2", -1 / 2)), sol)

        delayed_sin2(t) = t < 1 ? 0 : sin(2π * 1e-6 * Ω00 / 2 * (t - 1))^2
        ex_ionsim_tdE = real.(expect(ionprojector(T, ("D5/2", -1 / 2)), sol))
        ex_analyt_tdE = @.delayed_sin2(tout)
        @test isapprox(ex_ionsim_tdE, ex_analyt_tdE, rtol = 1e-5)

        # check that δν(t) shifts sideband frequencies appropriately
        L.E = E
        tspan = 0:1e-3:30
        mode.δν = t -> 20e3
        L.Δ = mode.ν + 20e3
        h = hamiltonian(T, timescale=1e-6, rwa_cutoff = 1e5)
        tout, sol = timeevolution.schroedinger_dynamic(
            tspan,
            ionstate(T, ("S1/2", -1 / 2)) ⊗ mode[1],
            h
        )
        η = lambdicke(mode, L, C)
        ex_ionsim_δν = real.(expect(ionprojector(T, ("D5/2", -1 / 2)), sol))
        ex_analyt_δν = @.(sin(2π * sqrt(2) * (η / sqrt(1.02)) * Ω00 / 2 * 1e-6 * tout)^2)
        @test isapprox(ex_ionsim_δν, ex_analyt_δν, rtol = 1e-1)

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
            com_frequencies = (x = 3e6, y = 3e6, z = 1e6),
            selected_modes = (; z = [1])
        )
        T = Chamber(iontrap = chain, B = 4e-4, Bhat = (x̂ + ẑ) / √2, lasers = [L1, L2])
        L1.λ = transitionwavelength(C, (("S1/2", -1 / 2), ("D5/2", -1 / 2)), T)
        L2.λ = transitionwavelength(C, (("S1/2", -1 / 2), ("D5/2", -1 / 2)), T)
        mode = T.iontrap.selected_modes.z[1]
        # Set the laser parameters
        ϵ = 40e3
        d = 80  # corrects for AC stark shift from single-photon coupling to sidebands
        L1.Δ = mode.ν + ϵ - d
        L1.k = ẑ
        L1.ϵ = x̂

        L2.Δ = -mode.ν - ϵ + d
        L2.k = ẑ
        L2.ϵ = x̂

        mode.N = 15
        η = abs(lambdicke(mode, L1, C))
        Ω = √(1e3 * ϵ) / η  # This will give a 1kHz MS strength, since coupling goes like (ηΩ)^2/ϵ

        Efield_from_rabi_frequency!(Ω, T, 1, 1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)))
        Efield_from_rabi_frequency!(Ω, T, 2, 1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)))
        ψi = ionstate(T, ("S1/2", -1 / 2), ("S1/2", -1 / 2)) ⊗ mode[0]  # initial state
        h = hamiltonian(T, timescale=1e-6, rwa_cutoff = 5e5)
        tspan = 0:0.25:1000
        tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h)
        SS = expect(ionprojector(T, ("S1/2", -1 / 2), ("S1/2", -1 / 2)), sol)
        DD = expect(ionprojector(T, ("D5/2", -1 / 2), ("D5/2", -1 / 2)), sol)
        ex = analytical.two_ion_ms(tspan, 1e-6Ω, 1e-6mode.ν, 1e-6mode.ν + 1e-6ϵ, η, 0)
        @test isapprox(ex[1], SS, rtol = 1e-2)
        @test isapprox(ex[2], DD, rtol = 1e-2)

        # add checks for fast and hot gates
    end
end  # end suppress
