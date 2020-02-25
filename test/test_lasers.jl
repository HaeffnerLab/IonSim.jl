using Test, IonSim

@testset "lasers" begin
    @test_throws AssertionError laser(ϵ=(x=2, y=0, z=0))
    @test_throws AssertionError laser(k=(x=2, y=0, z=0))
    @test_throws AssertionError laser(ϵ=(x=1, y=0, z=0), k=(x=1, y=0, z=0))
    @test_throws AssertionError laser(pointing=[(1, 0.5), (1, 0.5)])
    @test_throws AssertionError laser(pointing=[(1, 2.0)])
    t = 0:1:100
    L = laser(E=1, ϕ=1)
    ones = [1 for _ in t]
    @assert L.E.(t) == ones && L.ϕ.(t) == ones
    L = laser(E=sin, ϕ=sin)
    @assert L.E.(t) == sin.(t) && L.ϕ.(t) == sin.(t)
    @test_throws AssertionError L.ϵ = (x=2, y=0, z=0)
    @test_throws AssertionError L.k = (x=2, y=0, z=0)
    @test_throws AssertionError L.pointing = [(1, 0.5), (1, 0.5)]
    @test_throws AssertionError L.pointing = [(1, 2.0)]
    L.λ = 1
    @test L.λ == 1
end