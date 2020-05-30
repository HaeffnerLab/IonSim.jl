using QuantumOptics: NLevelBasis, CompositeBasis, FockBasis
using Test, IonSim


@testset "traps-required_fields" begin
    C = ca40()
    L1 = laser(); L1.pointing = [(1, 1.0), (2, 1.0)]
    L2 = laser(); L2.pointing = [(1, 1.0), (2, 1.0)]
    lasers = [L1, L2]
    chain = linearchain(ions=[C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[1], y=[], z=[1]))
    T = trap(label="mytrap", configuration=chain, B=6e-4, ∇B=2, Bhat=(x̂ + ẑ)/√2, δB=1, lasers=lasers);
    @test label(T) == "mytrap"
    @test configuration(T) == chain
    @test Bhat(T) == (x̂ + ẑ)/√2
    @test gradB(T) == 2
    @test get_lasers(T) == lasers
    t = 0:1:100
    ones = [1 for _ in t]
    @test deltaB(T).(t) == ones
    T.δB = sin
    @test deltaB(T).(t) == sin.(t)
    @test typeof(get_basis(T)) == CompositeBasis{Array{Int64,1},Tuple{NLevelBasis{Int64},NLevelBasis{Int64},FockBasis{Int64},FockBasis{Int64}}}
end


@testset "traps-main" begin
    C = ca40()
    L1 = laser(); L1.pointing = [(1, 1.0), (2, 1.0)]
    L2 = laser(); L2.pointing = [(1, 1.0), (2, 1.0)]
    lasers = [L1, L2]
    chain = linearchain(ions=[C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[1], y=[], z=[1]))
    T = trap(configuration=chain, lasers=[L1, L1]) 
    @test T.lasers[1] !== T.lasers[2]
    @test_throws AssertionError trap(configuration=chain, lasers=[L1], Bhat=(x=2, y=0, z=0))
    L3 = laser(); L3.pointing = [(3, 1.0)]
    @test_throws AssertionError trap(configuration=chain, lasers=[L3])
    C.label = "ion1"
    C2 = ca40(); C2.label = "ion2"
    chain.vibrational_modes.x[1].label = "radial"
    chain.vibrational_modes.z[1].label = "axial"
    T = trap(label="mytrap", configuration=chain, B=6e-4, ∇B=2, Bhat=(x̂ + ẑ)/√2, δB=1, lasers=lasers);
    @test T["ion1"] == C
    @test T["axial"] == chain.vibrational_modes.z[1]
end


@testset "traps-general" begin
    C = ca40()
    C1 = copy(C)
    L = laser(); L.pointing = [(1, 1.0)]
    lasers = [L]
    chain = linearchain(ions=[C, C1], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[], y=[], z=[1]))
    T = trap(configuration=chain, lasers=[L], Bhat=(x̂ + ŷ + ẑ)/√3) 
    
    # Efield_from_pi_time
    @test Efield_from_pi_time(1e-6, (x̂ + ŷ + ẑ)/√3, L, T.configuration.ions[1], ("S-1/2", "D-1/2")) ≈ 153739.56 rtol=1e-4
    @test Efield_from_pi_time(1e-6, (x̂ + ŷ + ẑ)/√3, L, T.configuration.ions[1], ("S-1/2", "D-1/2")) == Efield_from_pi_time(1e-6, T, 1, 1, ("S-1/2", "D-1/2"))
    Efield_from_pi_time!(1e-6, (x̂ + ŷ + ẑ)/√3, L, T.configuration.ions[1], ("S-1/2", "D-1/2"))
    @test T.lasers[1].E(0) == Efield_from_pi_time(1e-6, (x̂ + ŷ + ẑ)/√3, L, T.configuration.ions[1], ("S-1/2", "D-1/2"))

    # Efield_from_rabi_frequency
    @test T.lasers[1].E(0) ≈ Efield_from_rabi_frequency(1 / (2 * 1e-6), (x̂ + ŷ + ẑ)/√3, L, T.configuration.ions[1], ("S-1/2", "D-1/2"))
    Efield_from_rabi_frequency!(5e-3, (x̂ + ŷ + ẑ)/√3, L, T.configuration.ions[1], ("S-1/2", "D-1/2"))
    @test T.lasers[1].E(0) == Efield_from_rabi_frequency(5e-3, (x̂ + ŷ + ẑ)/√3, L, T.configuration.ions[1], ("S-1/2", "D-1/2"))

    # transition_frequency
    T.B = 4e-4
    @test transition_frequency(T, C, ("S-1/2", "D-1/2")) ≈ 2.2393e6 rtol=1e-4

    # set_gradient
    f = transition_frequency(T, C, ("S-1/2", "D-1/2"))
    set_gradient!(T, (C.number, C1.number), ("S-1/2", "D-1/2"), 1e6)
    @test transition_frequency(T, C, ("S-1/2", "D-1/2")) ≈ f - 1e6 / 2 rtol=1e-4
end
