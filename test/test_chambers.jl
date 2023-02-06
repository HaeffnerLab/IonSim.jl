using QuantumOptics: NLevelBasis, CompositeBasis, FockBasis
using Test, IonSim
using IonSim.PhysicalConstants: ħ, c
using Suppressor

@suppress_err begin

    # set up system
    C = Ca40()
    λ = transitionwavelength(C, ("S1/2", "D5/2"))
    L1 = Laser(λ=λ)
    L2 = Laser(λ=λ)
    chain = LinearChain(
        ions=[C, C],
        comfrequencies=(x=3e6, y=3e6, z=1e6),
        selectedmodes=(; z=[1])
    )

    @testset "chambers -- Chamber" begin
        T = Chamber(iontrap=chain, lasers=[L1, L1])

        # test construction of CompositeBasis
        @test basis(T).bases[1] ≡ T.iontrap.ions[1]
        @test basis(T).bases[2] ≡ T.iontrap.ions[2]
        @test basis(T).bases[3] ≡ T.iontrap.selectedmodes.z[1]

        # test construction of T.δB
        t = 0:1:100
        ones = [1 for _ in t]
        bfield_fluctuation!(T, 1)
        @test T.δB.(t) == ones
        @test T._cnst_δB
        bfield_fluctuation!(T, sin)
        @test T.δB.(t) == sin.(t)
        @test !T._cnst_δB

        # test for warning when lasers=[L, L, L, ...] where L point to the same thing
        warn = "Some lasers point to the same thing. Making copies."
        @test_logs (:warn, warn) Chamber(iontrap=chain, lasers=[L1, L1, L1])
        # and test that, in this case, copies are made
        T1 = Chamber(iontrap=chain, lasers=[L1, L1, L1, L1], δB=t -> t)
        for i in 1:4, j in (i+1):4
            @test !(T1.lasers[i] ≡ T1.lasers[j])
        end

        # test setproperty for Chamber
        # setting T.δB to a constant
        bfield_fluctuation!(T, 100)
        @test T.δB(0) == 100
        # # setting the :basis directly should have no effect
        # basis = basis(T)
        # basis(T) = 100
        # @test basis(T) == basis
        # changing the iontrap should also update the basis
        chain1 = LinearChain(
            ions=[C, C, C],
            comfrequencies=(x=3e6, y=3e6, z=1e6),
            selectedmodes=(; z=[1])
        )
        @test length(basis(T).bases) ≡ 3
        iontrap!(T, chain1)
        @test T.iontrap ≡ chain1
        @test length(basis(T).bases) ≡ 4

        # make sure print/show doesn't throw errors
        print(T)
        show(T)
    end

    @testset "chambers -- general functions" begin
        Bhat = (x̂ + ŷ + ẑ) / √3
        T = Chamber(iontrap=chain, lasers=[L1], Bhat=Bhat)

        # globalbeam!
        L = Laser(λ=λ)
        globalbeam!(L, T)
        @test L.pointing == [(1, 1.0), (2, 1.0)]

        # intensity_from_pitime
        ion_index = 1
        laser_index = 1
        ion = T.iontrap.ions[ion_index]
        transition = (("S1/2", -1 / 2), ("D5/2", -1 / 2))
        # compare to specific pre-computed value for both methods
        I1 = intensity_from_pitime(L1, 1e-6, ion, transition, Bhat)
        I2 = intensity_from_pitime(laser_index, 1e-6, ion_index, transition, T)
        @test I1 ≈ 1.8556916e7 rtol = 1e-2
        @test I1 == I2
        # confirm in-place versions work
        intensity_from_pitime!(L1, 1e-6, ion, transition, Bhat)
        @test L1.I(0) == I1
        intensity_from_pitime!(laser_index, 1e-6, ion_index, transition, T)
        @test L1.I(0) == I1
        # shouldn't be able to have a laser argument where laser.pointing = []
        L = Laser(λ=λ)
        @test_throws AssertionError intensity_from_pitime(L, 1e-6, ion, transition, Bhat)
        L.pointing = [(1, 1.0)]
        @test isinf(intensity_from_pitime(L, 1e-6, ion, transition, x̂))

        # intensity_from_rabifrequency
        ion_index = 1
        laser_index = 1
        ion = T.iontrap.ions[ion_index]
        transition = (("S1/2", -1 / 2), ("D5/2", -1 / 2))
        # compare to specific pre-computed value for both methods
        I1 = intensity_from_rabifrequency(L1, 5e5, ion, transition, Bhat)
        I2 = intensity_from_rabifrequency(laser_index, 5e5, ion_index, transition, T)
        @test I1 ≈ 1.8556916e7 rtol = 1e-2
        @test I1 == I2
        # confirm in-place versions work
        intensity_from_rabifrequency!(L1, 5e5, ion, transition, Bhat)
        @test L1.I(0) == I1
        intensity_from_rabifrequency!(laser_index, 5e5, ion_index, transition, T)
        @test L1.I(0) == I1
        # shouldn't be able to have a laser argument where laser.pointing = []
        L = Laser(λ=λ)
        @test_throws AssertionError intensity_from_rabifrequency(
            L,
            5e5,
            ion,
            transition,
            Bhat
        )

        # transitionfrequency (test against pre-computed values)
        T.B = 4e-4
        f = transitionfrequency(C, transition, T)
        @test f ≈ 4.111550340833542e14
        @test transitionfrequency(C, transition, B=T.B) == f
        @test transitionfrequency(1, transition, T) == f

        # set_gradient
        bgradient!(T, (1, 2), transition, 1e6)
        f1 = transitionfrequency(T.iontrap.ions[1], transition, T)
        f2 = transitionfrequency(T.iontrap.ions[2], transition, T)
        @test abs(f1 - f2) ≈ 1e6

        # test :(==)
        chain1 = LinearChain(
            ions=[C, C],
            comfrequencies=(x=3e6, y=3e6, z=1e6),
            selectedmodes=(x=[1], z=[1])
        )
        T = Chamber(iontrap=chain1, lasers=[L1, L2])
        xmode = chain1.selectedmodes.x[1]
        zmode = chain1.selectedmodes.z[1]
        cb = basis(T)
        @test cb == C ⊗ C ⊗ xmode ⊗ zmode
        @test cb != C ⊗ C ⊗ zmode ⊗ xmode
        @test cb != C ⊗ C ⊗ C ⊗ xmode ⊗ zmode

        # test lambdicke
        # η = |k|cos(θ) * √(ħ / (2M ⋅ N ⋅ 2πν)); cos(θ) ≡ k̂ ⋅ mode_axis; N ≡ number of ions
        η(ν) = (2π / λ) * sqrt(ħ / (2 * mass(C) * 2 * 2π * ν))
        @test abs(lambdicke(zmode, C, L)) ≈ η(zmode.ν)
        L.k = (x̂ + ẑ) / √2
        @test abs(lambdicke(xmode, C, L)) ≈ η(xmode.ν) / √2
    end
end  # end suppress
