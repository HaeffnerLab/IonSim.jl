using QuantumOptics: NLevelBasis, nlevelstate
using Test, IonSim


@testset "ions-required_fields" begin
    # ca40
    C = ca40()
    @test typeof(label(C)) == String
    @test typeof(mass(C)) <: Real
    @test typeof(level_structure(C)) == OrderedDict{String,NamedTuple}
    @test typeof(selected_level_structure(C)) == OrderedDict{String,NamedTuple}
    @test typeof(matrix_elements(C)) == OrderedDict{Tuple,Function}
    @test typeof(get_basis(C)) <: NLevelBasis
    @test typeof(ion_number(C)) <: Nothing
    @test typeof(ion_position(C)) <: Nothing
    @test typeof(stark_shift(C)) == OrderedDict{String,Real}
end

@testset "ions--general" begin
    @test gJ(2, 5/2) == 6/5
    @test zeeman_shift(1e-4, 2, 5/2, 5/2) == 4.198873488465285e6

    # ca40
    C = ca40()
    C.stark_shift["S-1/2"] = 10
    zero_stark_shift(C)
    @test sum(values(C.stark_shift)) == 0
    @test_throws AssertionError C[""]
    @test C["S-1/2"] == nlevelstate(C.basis, 1)
    @test_throws AssertionError C.selected_level_structure = []
end

@testset "ca40" begin
    levels = ["S-1/2", "S+1/2", "D-1/2"]
    C = ca40(selected_level_structure=levels)
    @test all(keys(C.selected_level_structure) .== levels)
    @test all(keys(C.selected_matrix_elements) .== [("S-1/2", "D-1/2"), ("S+1/2", "D-1/2")])
    @test C.selected_matrix_elements[("S-1/2", "D-1/2")](1, 0, 45) ≈ 4.2248 rtol=1e-4
    @test sum(C.selected_matrix_elements[("S-1/2", "D-1/2")](1, 90, 0)) == 0
end