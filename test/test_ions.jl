using QuantumOptics: NLevelBasis, nlevelstate
using Test, IonSim
using IonSim.PhysicalConstants
using Suppressor
using Unitful
using InteractiveUtils

@suppress_err begin
    @testset "ions -- general" begin
        C = Ca40()
        # test for required fields
        @test typeof(speciesproperties(C)) <: NamedTuple
        @test typeof(sublevels(C)) == Vector{Tuple{String, Real}}
        @test typeof(sublevel_aliases(C)) == Dict{String, Tuple}
        @test isempty(sublevel_aliases(C))
        @test typeof(shape(C)) == Vector{Int64}
        @test typeof(stark_shift(C)) == OrderedDict{Tuple, INVERSE_TIME}
        @test typeof(ionnumber(C)) == Missing
        @test typeof(ionposition(C)) == Missing

        # test ==
        C1 = Ca40()
        @test C1 == C

        # test aliases
        set_sublevel_alias!(C, ("S1/2", -1 / 2), "S")
        @test sublevel_aliases(C) == Dict("S" => ("S1/2", -1 / 2))
        clear_sublevel_alias!(C, ("S1/2", -1 / 2))

        set_sublevel_alias!(
            C,
            [
                (("S1/2", -1 / 2), "S")
                (("D5/2", -1 / 2), "D")
            ]
        )
        @test sublevel_aliases(C) == Dict("S" => ("S1/2", -1 / 2), "D" => ("D5/2", -1 / 2))
        clear_sublevel_alias!(C, ["S", "D"])

        set_sublevel_alias!(C, Dict("0" => ("S1/2", 1 / 2), "1" => ("D5/2", 5 / 2)))
        @test sublevel_aliases(C) == Dict("0" => ("S1/2", 1 / 2), "1" => ("D5/2", 5 / 2))
        clear_all_sublevel_aliases!(C)

        # set some aliases for convenience
        set_sublevel_alias!(C, ("S1/2", -1 / 2), "S")
        set_sublevel_alias!(C, ("D5/2", -5 / 2), "D")

        #test stark shift
        set_stark_shift!(C, "S", 10.0u"1/s")
        @test stark_shift(C, "S") == 10.0u"1/s"
        zero_stark_shift!(C)
        @test stark_shift(C, "S") == 0.0u"1/s"

        # test levels and sublevels
        @test levels(C) == ["S1/2", "D5/2"]
        @test sublevelalias(C, ("S1/2", -1 / 2)) == "S"
        @test sublevelalias(C, ("D5/2", -1 / 2)) == nothing
        @test alias2sublevel(C, "S") == ("S1/2", -1 / 2)
        @test sublevel2level(C, "S") == "S1/2"
        @test sublevel2level(C, ("D5/2", -3 / 2)) == "D5/2"

        # test level and sublevel properties
        @test quantumnumbers(C, "D").m == -5 // 2
        @test quantumnumbers(C, ("D5/2", 1 / 2)).m == 1 // 2

        @test energy(C, "S1/2") == energy(C, "S")
        @test energy(C, "S1/2") != energy(C, "S", B = 1e-4u"T")
        @test transitionwavelength(C, ("S", "D")) ≈ 729.147u"nm" |> u"m"
    end

    @testset "ions -- species" begin
        # attempt to instantiate all Ion subtypes (use default sublevel selection)
        species = subtypes(Ion)
        for s in species
            ion = s()
            @test typeof(ion) <: Ion
        end
    end

    @testset "ions -- Ca40" begin
        C = Ca40()

        # set some aliases for convenience
        set_sublevel_alias!(C, ("S1/2", -1 / 2), "S")
        set_sublevel_alias!(C, ("D5/2", -5 / 2), "D")

        # test for general species properties
        @test mass(C) ≈ 6.635943757345042e-26u"kg"
        @test charge(C) ≈ e
        @test nuclearspin(C) == 0

        @test landegf(quantumnumbers(C, "D")) == 6 // 5
        @test landegf(C, "D") == 6 // 5

        # test zeeman shift using quantum numbers of the S and D states
        @test zeeman_shift(1e-4u"T", quantumnumbers(C, "S")) ≈ -1.3996244961550953e6u"1/s"
        @test zeeman_shift(1e-4u"T", quantumnumbers(C, "D")) ≈ -4.198873488465285e6u"1/s"

        # test zeeman shift using the ion itself as an input, which will use the custom-defined g-factors
        # for the S1/2 and D5/2 states and thus give slightly different results
        @test zeeman_shift(C, "S", 1e-4u"T") ≈ -1.4012037204665968e6u"1/s"
        @test zeeman_shift(C, "D", 1e-4u"T") ≈ -4.200042174919575e6u"1/s"

        @test leveltransitions(C) == [("S1/2", "D5/2")]
        @test length(subleveltransitions(C)) == 10

        @test lifetime(C, "D5/2") ≈ 1.16795141322121u"s"

        @test matrix_element(C, ("S", "D"), 1e5u"V/m", x̂, ŷ, ẑ) ≈
              472761.18184781645u"1/s"

        # # make sure improper indexing of Ca40 yields an AssertionError
        @test_throws AssertionError C[""]

        # # test indexing
        C1 = Ca40([("S1/2", -1 / 2), ("D5/2", -5 / 2)])
        set_sublevel_alias!(C1, ("S1/2", -1 / 2), "S")
        set_sublevel_alias!(C1, ("D5/2", -5 / 2), "D")
        @test C1[("S1/2", -1 / 2)].data == ComplexF64[1; 0]
        @test C1[("D5/2", -5 / 2)].data == ComplexF64[0; 1]
        @test C1["S"].data == ComplexF64[1; 0]
        @test C1["D"].data == ComplexF64[0; 1]

        # # test set properties
        # @test_throws AssertionError C.selected_level_structure = []

        # # test get properties
        warning = "ion has not been added to a configuration"
        @test_logs (:warn, warning) C.ionnumber
        @test_logs (:warn, warning) C.position
    end
end  # end suppress
