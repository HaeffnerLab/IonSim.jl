using Test, IonSim, IonSim.PhysicalConstants, Suppressor, Unitful

@testset "constants -- PhysicalConstant" begin
    # test that algebraic operations on physical constants return a unitful quantity
    @test typeof(c - c) <: typeof(u"c0")
    @test typeof(c - 1u"m/s") <:  typeof(u"c0")
    @test typeof(1u"m/s" - c) <:  typeof(u"c0")
    @test typeof(c + c) <:  typeof(u"c0")
    @test typeof(1u"m/s" + c) <:  typeof(u"c0")
    @test typeof(c + 1u"m/s") <:  typeof(u"c0")
    @test typeof(c * c) <:  typeof(u"c0"^2)
    @test typeof(1 * c) <: typeof(u"c0")
    @test typeof(c * 1) <: typeof(u"c0")
    @test typeof(c / c) <: Number
    @test typeof(1 / c) <: typeof(u"c0"^-1)
    @test typeof(c / 1) <: typeof(u"c0"*1.0) # this becomes a float
    @test typeof(c^3) <: typeof(u"c0"^3)
    @test typeof(2^α) <: Number
    @test typeof(sqrt(c)) <: typeof(u"c0"^0.5)
    @test typeof(α^α) <: Number

    # Test comparisons for physical constants
    @test ustrip(c) != c
    @test c != ustrip(c)
    @test c == c
    @test c > c - 1u"m/s"
    @test c > c - 1u"m/s"
    @test ustrip(c) < ustrip(c) + 1
    @test c < c + 1u"m/s"
    @test c ≥ c
    @test ustrip(c) ≤ ustrip(c)
    @test c ≤ c
    @test α < ustrip(c)
    @test α > ustrip(ħ)

    # test print/show for physicalconstants
    out = @capture_out print(c)
    @test out == "$(ustrip(c)) $(unit(c))"

    @test IonSim._print_axis(x̂) == "x̂"
    @test IonSim._print_axis(ŷ) == "ŷ"
    @test IonSim._print_axis(ẑ) == "ẑ"
    @test IonSim._print_axis((x = 1, y = 1, z = 1)) == string((x = 1, y = 1, z = 1))
end

@testset "constants -- other constants" begin
    @test x̂ / 2 == (x = 1 / 2, y = 0, z = 0)
    @test 2x̂ == (x = 2, y = 0, z = 0)
    @test x̂ * 2 == (x = 2, y = 0, z = 0)
    @test x̂ + x̂ == (x = 2, y = 0, z = 0)
    @test x̂ - x̂ == (x = 0, y = 0, z = 0)
    @test ndot(x̂, x̂) == 1
    @test ndot(x̂, ŷ) == 0
end
