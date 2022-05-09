using QuantumOptics: NLevelBasis, CompositeBasis, FockBasis
using Test, IonSim
using IonSim.PhysicalConstants: ħ, c
using Suppressor
using Unitful

@suppress_err begin

    # set up system
    C = Ca40()
    λ = transitionwavelength(C, ("S1/2", "D5/2"))
    L1 = Laser(λ = λ)
    L2 = Laser(λ = λ)
    chain = LinearChain(
        ions = [C, C],
        com_frequencies = (x = 3e6u"1/s", y = 3e6u"1/s", z = 1e6u"1/s"),
        vibrational_modes = (; z = [1])
    )

    @testset "traps -- Trap" begin
        T = Trap(configuration = chain, lasers = [L1, L1])

        # test construction of CompositeBasis
        @test T.basis.bases[1] ≡ T.configuration.ions[1]
        @test T.basis.bases[2] ≡ T.configuration.ions[2]
        @test T.basis.bases[3] ≡ T.configuration.vibrational_modes.z[1]

        # test construction of T.δB
        t = 0:1:100
        ones = [1 for _ in t]
        T.δB = 1u"T"
        @test T.δB.(t) == 1u"T" * ones
        @test T._cnst_δB
        T.δB = x -> 1u"T" * sin.(x) #TODO: do I need to have time be marked?
        @test T.δB.(t) == 1u"T" * sin.(t)
        @test !T._cnst_δB

        # test for warning when lasers=[L, L, L, ...] where L point to the same thing
        warn = "Some lasers point to the same thing. Making copies."
        @test_logs (:warn, warn) Trap(configuration = chain, lasers = [L1, L1, L1])
        # and test that, in this case, copies are made
        T1 = Trap(configuration = chain, lasers = [L1, L1, L1, L1], δB = t -> t)
        for i in 1:4, j in (i + 1):4
            @test !(T1.lasers[i] ≡ T1.lasers[j])
        end

        # test setproperty for Trap
        # setting T.δB to a constant
        T.δB = 100u"T"
        @test T.δB(0) == 100u"T"
        # setting the :basis directly should have no effect
        basis = T.basis
        T.basis = 100
        @test T.basis == basis
        # changing the configuration should also update the basis
        chain1 = LinearChain(
            ions = [C, C, C],
            com_frequencies = (x = 3e6u"1/s", y = 3e6u"1/s", z = 1e6u"1/s"),
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

    @testset "traps -- general functions" begin
        Bhat = (x̂ + ŷ + ẑ) / √3
        T = Trap(configuration = chain, lasers = [L1], Bhat = Bhat)

        # get_basis
        @test T.basis == get_basis(T)

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
        E1 = Efield_from_pi_time(1e-6u"s", Bhat, L1, ion, transition)
        E2 = Efield_from_pi_time(1e-6u"s", T, laser_index, ion_index, transition)
        @test E1 ≈ 118245.11u"V/m" rtol = 1e-2
        @test E1 == E2
        # confirm in-place versions work
        Efield_from_pi_time!(1e-6u"s", Bhat, L1, ion, transition)
        @test L1.E(0) == E1
        Efield_from_pi_time!(1e-6u"s", T, laser_index, ion_index, transition)
        @test L1.E(0) == E1
        # shouldn't be able to have a laser argument where laser.pointing = []
        L = Laser(λ = λ)
        @test_throws AssertionError Efield_from_pi_time(1e-6u"s", Bhat, L, ion, transition)
        L.pointing = [(1, 1.0)]
        @test isinf(Efield_from_pi_time(1e-6u"s", x̂, L, ion, transition))

        # Efield_from_rabi_frequency
        ion_index = 1
        laser_index = 1
        ion = T.configuration.ions[ion_index]
        transition = (("S1/2", -1 / 2), ("D5/2", -1 / 2))
        # compare to specific pre-computed value for both methods
        E1 = Efield_from_rabi_frequency(5e5u"1/s", Bhat, L1, ion, transition)
        E2 = Efield_from_rabi_frequency(5e5u"1/s", T, laser_index, ion_index, transition)
        @test E1 ≈ 118245.11u"V/m" rtol = 1e-2
        @test E1 == E2
        # confirm in-place versions work
        Efield_from_rabi_frequency!(5e5u"1/s", Bhat, L1, ion, transition)
        @test L1.E(0) == E1
        Efield_from_rabi_frequency!(5e5u"1/s", T, laser_index, ion_index, transition)
        @test L1.E(0) == E1
        # shouldn't be able to have a laser argument where laser.pointing = []
        L = Laser(λ = λ)
        @test_throws AssertionError Efield_from_rabi_frequency(
            5e5u"1/s",
            Bhat,
            L,
            ion,
            transition
        )

        # transitionfrequency (test against pre-computed values)
        T.B = 4e-4u"T"
        f = transitionfrequency(C, transition, T)
        @test f ≈ 4.111550340833542e14u"1/s"
        @test transitionfrequency(C, transition, B = T.B) == f
        @test transitionfrequency(1, transition, T) == f

        # set_gradient
        set_gradient!(T, (1, 2), transition, 1e6u"1/s")
        f1 = transitionfrequency(T.configuration.ions[1], transition, T)
        f2 = transitionfrequency(T.configuration.ions[2], transition, T)
        @test abs(f1 - f2) ≈ 1e6u"1/s"

        # test :(==)
        chain1 = LinearChain(
            ions = [C, C],
            com_frequencies = (x = 3e6u"1/s", y = 3e6u"1/s", z = 1e6u"1/s"),
            vibrational_modes = (x = [1], z = [1])
        )
        T = Trap(configuration = chain1, lasers = [L1, L2])
        xmode = chain1.vibrational_modes.x[1]
        zmode = chain1.vibrational_modes.z[1]
        cb = get_basis(T)
        @test cb == C ⊗ C ⊗ xmode ⊗ zmode
        @test cb != C ⊗ C ⊗ zmode ⊗ xmode
        @test cb != C ⊗ C ⊗ C ⊗ xmode ⊗ zmode

        # test get_η
        # η = |k|cos(θ) * √(ħ / (2M ⋅ N ⋅ 2πν)); cos(θ) ≡ k̂ ⋅ mode_axis; N ≡ number of ions
        η(ν) = (2π / λ) * sqrt(ħ / (2 * mass(C) * 2 * 2π * ν)) |> NoUnits
        @test abs(get_η(zmode, L, C)) ≈ η(zmode.ν)
        L.k = (x̂ + ẑ) / √2
        @test abs(get_η(xmode, L, C)) ≈ η(xmode.ν) / √2
        @test typeof(η(xmode.ν)) <: Number
    end
end  # end suppress
