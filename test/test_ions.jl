using QuantumOptics: NLevelBasis, nlevelstate
using Test, IonSim
using IonSim.Properties: IonProperties
using IonSim.PhysicalConstants
using Suppressor

@suppress_err begin
    @testset "ions -- general" begin
        C = Ion(CA40_PROPERTIES, nothing)
        # test for required fields
        @test typeof(speciesproperties(C)) <: IonProperties
        @test typeof(sublevels(C)) == Vector{Tuple{String, Real}}
        @test typeof(sublevelaliases(C)) == Dict{String, Tuple}
        @test isempty(sublevelaliases(C))
        @test typeof(shape(C)) == Vector{Int64}
        @test typeof(manualshift(C)) == OrderedDict{Tuple, Real}
        @test typeof(ionnumber(C)) == Missing
        @test typeof(ionposition(C)) == Missing

        # test ==
        C1 = Ion(CA40_PROPERTIES, nothing)
        @test C1 == C

        # test aliases
        sublevelalias!(C, ("S1/2", -1 / 2), "S")
        @test sublevelaliases(C) == Dict("S" => ("S1/2", -1 / 2))
        clearsublevelalias!(C, ("S1/2", -1 / 2))

        sublevelalias!(
            C,
            [
                (("S1/2", -1 / 2), "S")
                (("D5/2", -1 / 2), "D")
            ]
        )
        @test sublevelaliases(C) == Dict("S" => ("S1/2", -1 / 2), "D" => ("D5/2", -1 / 2))
        clearsublevelalias!(C, ["S", "D"])

        sublevelalias!(C, Dict("0" => ("S1/2", 1 / 2), "1" => ("D5/2", 5 / 2)))
        @test sublevelaliases(C) == Dict("0" => ("S1/2", 1 / 2), "1" => ("D5/2", 5 / 2))
        clearsublevelalias!(C)

        # set some aliases for convenience
        sublevelalias!(C, ("S1/2", -1 / 2), "S")
        sublevelalias!(C, ("D5/2", -5 / 2), "D")

        #test manual shift
        manualshift!(C, "S", 10.0)
        @test manualshift(C, "S") == 10.0
        zeromanualshift!(C)
        @test manualshift(C, "S") == 0.0

        # test levels and sublevels
        @test levels(C) == ["S1/2", "D5/2"]
        @test sublevelalias(C, ("S1/2", -1 / 2)) == "S"
        @test sublevelalias(C, ("D5/2", -1 / 2)) == nothing
        @test sublevel(C, "S") == ("S1/2", -1 / 2)
        @test level(C, "S") == "S1/2"
        @test level(C, ("D5/2", -3 / 2)) == "D5/2"

        # test level and sublevel properties
        @test quantumnumbers(C, "D").m == -5 // 2
        @test quantumnumbers(C, ("D5/2", 1 / 2)).m == 1 // 2

        @test energy(C, "S1/2") == energy(C, "S")
        @test energy(C, "S1/2") != energy(C, "S", B=1e-4)
        @test transitionwavelength(C, ("S", "D")) ≈ 7.29147e-7
    end

    @testset "ions -- Ca40" begin
        C = Ion(CA40_PROPERTIES, nothing)

        # set some aliases for convenience
        sublevelalias!(C, ("S1/2", -1 / 2), "S")
        sublevelalias!(C, ("D5/2", -5 / 2), "D")

        # test for general species properties
        @test mass(C) ≈ 6.635943757345042e-26
        @test charge(C) ≈ e
        @test nuclearspin(C) == 0

        @test landegf(quantumnumbers(C, "D")) == 6 // 5
        @test landegf(C, "D") == 6 // 5

        # test zeeman shift using quantum numbers of the S and D states
        @test zeemanshift(1e-4, quantumnumbers(C, "S")) ≈ -1.3996244961550953e6
        @test zeemanshift(1e-4, quantumnumbers(C, "D")) ≈ -4.198873488465285e6

        # test zeeman shift using the ion itself as an input, which will use the custom-defined g-factors
        # for the S1/2 and D5/2 states and thus give slightly different results
        @test zeemanshift(C, "S", 1e-4) ≈ -1.4012037204665968e6
        @test zeemanshift(C, "D", 1e-4) ≈ -4.200042174919575e6

        @test leveltransitions(C) == [("S1/2", "D5/2")]
        @test length(subleveltransitions(C)) == 10

        @test lifetime(C, "D5/2") ≈ 1.16795141322121

        @test matrixelement(C, ("S", "D"), 1.327209365e7, ŷ, x̂, ẑ) ≈ 472761.18184781645

        # make sure improper indexing of Ca40 yields an AssertionError
        @test_throws AssertionError C[""]

        # # test indexing
        C1 = Ion(CA40_PROPERTIES, [("S1/2", -1 / 2), ("D5/2", -5 / 2)])
        sublevelalias!(C1, ("S1/2", -1 / 2), "S")
        sublevelalias!(C1, ("D5/2", -5 / 2), "D")
        @test C1[("S1/2", -1 / 2)].data == ComplexF64[1; 0]
        @test C1[("D5/2", -5 / 2)].data == ComplexF64[0; 1]
        @test C1["S"].data == ComplexF64[1; 0]
        @test C1["D"].data == ComplexF64[0; 1]

        # test set properties
        # @test_throws AssertionError C.selected_level_structure = []

        # # test get properties
        warning = "ion has not been added to an IonTrap"
        @test_logs (:warn, warning) ionnumber(C)
        @test_logs (:warn, warning) ionposition(C)
    end
end  # end suppress
