using Test, IonSim, IonSim.PhysicalConstants, Suppressor

@testset "constants -- PhysicalConstant" begin
    # test that algebraic operations on physical constants return a <:Number
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
    @test typeof(sqrt(c)) <: Number
    @test typeof(α^α) <: Number

    # Test comparisons for physical constants
    @test c.x == c
    @test c == c.x
    @test c == c
    @test c.x > c - 1
    @test c > c.x - 1
    @test c > c - 1
    @test c.x < c + 1
    @test c < c.x + 1
    @test c < c + 1
    @test c.x ≥ c
    @test c ≥ c.x
    @test c ≥ c
    @test c.x ≤ c
    @test c ≤ c.x
    @test c ≤ c
    @test α < c
    @test α > ħ

    # test print/show for physicalconstants
    out = @capture_out print(c)
    @test out == "$(c.x) [$(c.units)]"

    @test IonSim._print_axis(x̂) == "x̂"
    @test IonSim._print_axis(ŷ) == "ŷ"
    @test IonSim._print_axis(ẑ) == "ẑ"
    @test IonSim._print_axis((x=1, y=1, z=1)) == string((x=1, y=1, z=1))
end

@testset "constants -- other constants" begin
    @test x̂ / 2 == (x=1 / 2, y=0, z=0)
    @test 2x̂ == (x=2, y=0, z=0)
    @test x̂ * 2 == (x=2, y=0, z=0)
    @test x̂ + x̂ == (x=2, y=0, z=0)
    @test x̂ - x̂ == (x=0, y=0, z=0)
    @test ndot(x̂, x̂) == 1
    @test ndot(x̂, ŷ) == 0
end
