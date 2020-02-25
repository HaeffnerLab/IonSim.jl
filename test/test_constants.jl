using Test, IonSim, IonSim.PhysicalConstants


@testset "physicalconstants" begin
    @test typeof(c - c) <: Number
    @test typeof(c - 1) <: Number
    @test typeof(1 - c) <: Number
    @test typeof(c + c) <: Number
    @test typeof(1 + c) <: Number
    @test typeof(c + 1) <: Number
    @test typeof(c * c) <: Number
    @test typeof(1 * c) <: Number
    @test typeof(c * 1) <: Number
    @test typeof(c / c) <: Number
    @test typeof(1 / c) <: Number
    @test typeof(c / 1) <: Number
    @test typeof(c^3) <: Number
    @test typeof(c^α) <: Number
    @test typeof(2^α) <: Number
    @test c.x == c
    @test c == c.x
    @test typeof(sqrt(c)) <: Number
end

@testset "other constants" begin
    @test x̂/2 == (x=1/2, y=0, z=0)
    @test 2x̂ == (x=2, y=0, z=0)
    @test x̂*2 == (x=2, y=0, z=0)
    @test x̂ + x̂ == (x=2, y=0, z=0)
    @test x̂ - x̂ == (x=0, y=0, z=0)
    @test ndot(x̂, x̂) == 1
    @test ndot(x̂, ŷ) == 0
end
