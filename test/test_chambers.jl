using QuantumOptics: NLevelBasis, CompositeBasis, FockBasis
using Test, IonSim
using IonSim.PhysicalConstants: ħ, c
using Suppressor

@suppress_err begin

    # set up system
    C = Ca40()
    λ = transitionwavelength(C, ("S1/2", "D5/2"))
    L1 = Laser(λ = λ)
    L2 = Laser(λ = λ)
    chain = LinearChain(
        ions = [C, C],
        com_frequencies = (x = 3e6, y = 3e6, z = 1e6),
        vibrational_modes = (; z = [1])
    )

    @testset "chambers -- Chamber" begin
        T = Chamber(configuration = chain, lasers = [L1, L1])

        # test construction of CompositeBasis
        @test T.basis.bases[1] ≡ T.configuration.ions[1]
        @test T.basis.bases[2] ≡ T.configuration.ions[2]
        @test T.basis.bases[3] ≡ T.configuration.vibrational_modes.z[1]

        # test construction of T.δB
        t = 0:1:100
        ones = [1 for _ in t]
        T.δB = 1
        @test T.δB.(t) == ones
        @test T._cnst_δB
        T.δB = sin
        @test T.δB.(t) == sin.(t)
        @test !T._cnst_δB

        # test for warning when lasers=[L, L, L, ...] where L point to the same thing
        warn = "Some lasers point to the same thing. Making copies."
        @test_logs (:warn, warn) Chamber(configuration = chain, lasers = [L1, L1, L1])
        # and test that, in this case, copies are made
        T1 = Chamber(configuration = chain, lasers = [L1, L1, L1, L1], δB = t -> t)
        for i in 1:4, j in (i + 1):4
            @test !(T1.lasers[i] ≡ T1.lasers[j])
        end

        # test setproperty for Chamber
        # setting T.δB to a constant
        T.δB = 100
        @test T.δB(0) == 100
        # setting the :basis directly should have no effect
        basis = T.basis
        T.basis = 100
        @test T.basis == basis
        # changing the configuration should also update the basis
        chain1 = LinearChain(
            ions = [C, C, C],
            com_frequencies = (x = 3e6, y = 3e6, z = 1e6),
            vibrational_modes = (; z = [1])
        )
        @test length(T.basis.bases) ≡ 3
        T.configuration = chain1
        @test T.configuration ≡ chain1
        @test length(T.basis.bases) ≡ 4

        # make sure print/show doesn't throw errors
        print(T)
        show(T)
    end

    @testset "chambers -- general functions" begin
        Bhat = (x̂ + ŷ + ẑ) / √3
        T = Chamber(configuration = chain, lasers = [L1], Bhat = Bhat)

        # basis
        @test T.basis == basis(T)

        # global_beam!
        L = Laser(λ = λ)
        global_beam!(T, L)
        @test L.pointing == [(1, 1.0), (2, 1.0)]

        # Efield_from_pi_time
        ion_index = 1
        laser_index = 1
        ion = T.configuration.ions[ion_index]
        transition = (("S1/2", -1 / 2), ("D5/2", -1 / 2))
        # compare to specific pre-computed value for both methods
        E1 = Efield_from_pi_time(1e-6, Bhat, L1, ion, transition)
        E2 = Efield_from_pi_time(1e-6, T, laser_index, ion_index, transition)
        @test E1 ≈ 118245.11 rtol = 1e-2
        @test E1 == E2
        # confirm in-place versions work
        Efield_from_pi_time!(1e-6, Bhat, L1, ion, transition)
        @test L1.E(0) == E1
        Efield_from_pi_time!(1e-6, T, laser_index, ion_index, transition)
        @test L1.E(0) == E1
        # shouldn't be able to have a laser argument where laser.pointing = []
        L = Laser(λ = λ)
        @test_throws AssertionError Efield_from_pi_time(1e-6, Bhat, L, ion, transition)
        L.pointing = [(1, 1.0)]
        @test isinf(Efield_from_pi_time(1e-6, x̂, L, ion, transition))

        # Efield_from_rabi_frequency
        ion_index = 1
        laser_index = 1
        ion = T.configuration.ions[ion_index]
        transition = (("S1/2", -1 / 2), ("D5/2", -1 / 2))
        # compare to specific pre-computed value for both methods
        E1 = Efield_from_rabi_frequency(5e5, Bhat, L1, ion, transition)
        E2 = Efield_from_rabi_frequency(5e5, T, laser_index, ion_index, transition)
        @test E1 ≈ 118245.11 rtol = 1e-2
        @test E1 == E2
        # confirm in-place versions work
        Efield_from_rabi_frequency!(5e5, Bhat, L1, ion, transition)
        @test L1.E(0) == E1
        Efield_from_rabi_frequency!(5e5, T, laser_index, ion_index, transition)
        @test L1.E(0) == E1
        # shouldn't be able to have a laser argument where laser.pointing = []
        L = Laser(λ = λ)
        @test_throws AssertionError Efield_from_rabi_frequency(
            5e5,
            Bhat,
            L,
            ion,
            transition
        )

        # transitionfrequency (test against pre-computed values)
        T.B = 4e-4
        f = transitionfrequency(C, transition, T)
        @test f ≈ 4.111550340833542e14
        @test transitionfrequency(C, transition, B = T.B) == f
        @test transitionfrequency(1, transition, T) == f

        # set_gradient
        set_gradient!(T, (1, 2), transition, 1e6)
        f1 = transitionfrequency(T.configuration.ions[1], transition, T)
        f2 = transitionfrequency(T.configuration.ions[2], transition, T)
        @test abs(f1 - f2) ≈ 1e6

        # test :(==)
        chain1 = LinearChain(
            ions = [C, C],
            com_frequencies = (x = 3e6, y = 3e6, z = 1e6),
            vibrational_modes = (x = [1], z = [1])
        )
        T = Chamber(configuration = chain1, lasers = [L1, L2])
        xmode = chain1.vibrational_modes.x[1]
        zmode = chain1.vibrational_modes.z[1]
        cb = basis(T)
        @test cb == C ⊗ C ⊗ xmode ⊗ zmode
        @test cb != C ⊗ C ⊗ zmode ⊗ xmode
        @test cb != C ⊗ C ⊗ C ⊗ xmode ⊗ zmode

        # test lambdicke
        # η = |k|cos(θ) * √(ħ / (2M ⋅ N ⋅ 2πν)); cos(θ) ≡ k̂ ⋅ mode_axis; N ≡ number of ions
        η(ν) = (2π / λ) * sqrt(ħ / (2 * mass(C) * 2 * 2π * ν))
        @test abs(lambdicke(zmode, L, C)) ≈ η(zmode.ν)
        L.k = (x̂ + ẑ) / √2
        @test abs(lambdicke(xmode, L, C)) ≈ η(xmode.ν) / √2
    end
end  # end suppress
