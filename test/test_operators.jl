using QuantumOptics: ⊗, number, expect
using Test, IonSim
using IonSim: thermalstate


@testset "operators" begin
    C = ca40()
    chain = linearchain(
        ions=[C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[1], y=[], z=[1])
    )
    mode = chain.vibrational_modes.z[1]
    mode.N = 100
    T = trap(configuration=chain)
    @test sigma(C, "S-1/2", "D-5/2") == C["S-1/2"] ⊗ C["D-5/2"]'
    @test ion_state(T, "S-1/2", "D-5/2") == C["S-1/2"] ⊗ C["D-5/2"]
    @test_throws AssertionError ion_state(T, "S-1/2")
    @test thermalstate(mode.basis, 5) == thermalstate(mode, 5)
    @test expect(number(mode.basis), thermalstate(mode, 5)) ≈ 5 rtol=1e-3
end