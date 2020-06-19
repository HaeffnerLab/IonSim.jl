using QuantumOptics: NLevelBasis, CompositeBasis, FockBasis
using Test, IonSim


# setup system
C = Ca40()
L1 = Laser(); L2 = Laser()
chain = LinearChain(
            ions=[C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), vibrational_modes=(;z=[1])
        )


@testset "traps -- Trap" begin
    T = Trap(configuration=chain, lasers=[L1, L1]) 
    
    # test construction of CompositeBasis
    @test T.basis.bases[1] ≡ T.configuration.ions[1]
    @test T.basis.bases[2] ≡ T.configuration.ions[2]
    @test T.basis.bases[3] ≡ T.configuration.vibrational_modes.z[1]

    # test construction of T.δB
    t = 0:1:100
    ones = [1 for _ in t]
    T.δB = 1
    @test T.δB.(t) == ones
    T.δB = sin
    @test T.δB.(t) == sin.(t)

    # test for warning when lasers=[L, L, L, ...] where L point to the same thing
    warn = "Some lasers point to the same thing. Making copies."
    @test_logs (:warn, warn) Trap(configuration=chain, lasers=[L1, L1, L1])
    # and test that, in this case, copies are made
    T1 = Trap(configuration=chain, lasers=[L1, L1, L1, L1])
    for i in 1:4, j in i+1:4
        @test !(T1.lasers[i] ≡ T1.lasers[j])
    end

    # test setproperty for Trap
    # setting T.δB to a constant
    T.δB = 100
    @test T.δB(0) == 100
    # setting the :basis directly should have no effect
    basis = T.basis
    T.basis = 100
    @test T.basis == basis
    # changing the configuration should also update the basis
    chain1 = LinearChain(
            ions=[C, C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), vibrational_modes=(;z=[1])
        )
    @test length(T.basis.bases) ≡ 3
    T.configuration = chain1
    @test T.configuration ≡ chain1
    @test length(T.basis.bases) ≡ 4
    
    # make sure print/show doesn't throw errors
    print(T); show(T)   
end


@testset "traps -- general functions" begin
    Bhat = (x̂ + ŷ + ẑ)/√3
    T = Trap(configuration=chain, lasers=[L1], Bhat=Bhat)

    # get_basis
    @test T.basis == get_basis(T)

    # global_beam!
    L = Laser()
    global_beam!(T, L)
    @test L.pointing == [(1, 1.0), (2, 1.0)]
    
    # Efield_from_pi_time
    ion_index = 1
    laser_index = 1
    ion = T.configuration.ions[ion_index]
    transition = ("S-1/2", "D-1/2")
    # compare to specific pre-computed value for both methods
    E1 = Efield_from_pi_time(1e-6, Bhat, L1, ion, transition)
    E2 = Efield_from_pi_time(1e-6, T, laser_index, ion_index, transition)
    @test E1 ≈ 153739.56 rtol=1e-2
    @test E1 == E2
    # confirm in-place versions work
    Efield_from_pi_time!(1e-6, Bhat, L1, ion, transition)
    @test L1.E(0) == E1
    Efield_from_pi_time!(1e-6, T, laser_index, ion_index, transition)
    @test L1.E(0) == E1
    # shouldn't be able to have a laser argument where laser.pointing = []
    L = Laser()
    @test_throws AssertionError Efield_from_pi_time(1e-6, Bhat, L, ion, transition)

    
    # Efield_from_rabi_frequency
    ion_index = 1
    laser_index = 1
    ion = T.configuration.ions[ion_index]
    transition = ("S-1/2", "D-1/2")
    # compare to specific pre-computed value for both methods
    E1 = Efield_from_rabi_frequency(5e5, Bhat, L1, ion, transition)
    E2 = Efield_from_rabi_frequency(5e5, T, laser_index, ion_index, transition)
    @test E1 ≈ 153739.56 rtol=1e-2
    @test E1 == E2
    # confirm in-place versions work
    Efield_from_rabi_frequency!(5e5, Bhat, L1, ion, transition)
    @test L1.E(0) == E1
    Efield_from_rabi_frequency!(5e5, T, laser_index, ion_index, transition)
    @test L1.E(0) == E1
    # shouldn't be able to have a laser argument where laser.pointing = []
    L = Laser()
    @test_throws AssertionError Efield_from_rabi_frequency(5e5, Bhat, L, ion, transition)

    # transition_frequency (test against pre-computed values)
    T.B = 4e-4
    f = transition_frequency(T, C, transition)
    @test f ≈ 2.2393e6 rtol=1e-4
    @test transition_frequency(T.B, C, transition) == f
    @test transition_frequency(T, 1, transition) == f

    # set_gradient
    set_gradient!(T, (1, 2), transition, 1e6)
    f1 = transition_frequency(T, T.configuration.ions[1], transition)
    f2 = transition_frequency(T, T.configuration.ions[2], transition)
    @test abs(f1 - f2) ≈ 1e6

    # test :(==)
    chain1 = LinearChain(
            ions=[C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), vibrational_modes=(x=[1], z=[1])
        )
    T = Trap(configuration=chain1, lasers=[L1, L2])
    xmode = chain1.vibrational_modes.x[1]
    zmode = chain1.vibrational_modes.z[1]
    cb = get_basis(T)
    @test cb == C ⊗ C ⊗ xmode ⊗ zmode
    @test cb != C ⊗ C ⊗ zmode ⊗ xmode
end
