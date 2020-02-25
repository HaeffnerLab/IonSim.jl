using Test, IonSim

@testset "vibrational_mode" begin
    t = 0:1:100
    ones = [1 for _ in t]
    vm = vibrational_mode("", 1e6, [1, 1], δν=1)
    @test vm.δν.(t) == ones

    vm.δν = sin
    @test vm.δν.(t) == sin.(t)

    vm.N = 20
    @test vm.basis.N == 19
    @test_throws AssertionError vm.N = 1.0
    @test_throws AssertionError vm.N = -1
end