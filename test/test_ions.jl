using QuantumOptics: NLevelBasis, nlevelstate
using Test, IonSim
const pc = IonSim.PhysicalConstants


@testset "Ca40" begin
    C = Ca40()
    # test for required fields
    @test mass(C) == pc.m_ca40
    @test typeof(level_structure(C)) == OrderedDict{String,NamedTuple}
    @test selected_level_structure(C) == level_structure(C)
    @test typeof(matrix_elements(C)) == OrderedDict{Tuple,Function}
    @test selected_matrix_elements(C) == matrix_elements(C)
    @test typeof(ion_number(C)) <: Nothing
    @test typeof(ion_position(C)) <: Nothing
    @test typeof(stark_shift(C)) == OrderedDict{String,Real}

    # test inner constructors for Ca40
    C1 = Ca40(["S-1/2", "D-1/2"])
    @test collect(keys(C1.selected_level_structure)) == ["S-1/2", "D-1/2"]
    C1copy = copy(C1)
    @test selected_level_structure(C1copy) == selected_level_structure(C1)

    # make sure print/show don't throw any errors
    print(C); show(C)

    # test matrix_element
    E = 100randn(); γ = 100randn(); ϕ = 100randn()
    @test matrix_element(C, ["S-1/2", "D-1/2"], E, γ, ϕ) == C.matrix_elements[("S-1/2", "D-1/2")](E, γ, ϕ)

    # Test matrix_element function for a specific collecion of parameters that have been
    # pre-evaluated independently. Hopefully this will fail if the function is constructed
    # improperly
    @test C.selected_matrix_elements[("S-1/2", "D-1/2")](1, 0, 45) ≈ 4.2248 rtol=1e-4
    @test sum(C.selected_matrix_elements[("S-1/2", "D-1/2")].(0, 90, 0:0.1:90)) == 0 
    
    # Test zero_stark_shift
    C.stark_shift["S-1/2"] = 10
    zero_stark_shift(C)
    @test sum(values(C.stark_shift)) == 0

    # make sure improper indexing of Ca40 yields an AssertionError
    @test_throws AssertionError C[""]

    # test indexing
    @test C1["S-1/2"].data == ComplexF64[1; 0]
    @test C1["D-1/2"].data == ComplexF64[0; 1]

    # test set properties
    @test_throws AssertionError C.selected_level_structure = []
    C.mass = 0  # shouldn't be allowed
    @test C.mass == pc.m_ca40
    C.selected_level_structure = ["S-1/2", "D-1/2"]
    @test C.selected_level_structure == C1.selected_level_structure

    # test get properties
    warning = "ion has not been added to a configuration"
    @test_logs (:warn, warning) C.number
    @test_logs (:warn, warning) C.position
end

@testset "ions -- general" begin

    # test specific case for Lande g-factor
    @test gJ(2, 5/2) == 6/5

    # test specific pre-evaluated case for zeeman_shift
    val = 4.198873488465285e6
    @test zeeman_shift(1e-4, 2, 5/2, 5/2) ≈ val
    @test zeeman_shift(1e-4, (2, 5/2, 5/2)) ≈ val
    @test zeeman_shift(B=1e-4, l=2, j=5/2, mⱼ=5/2) ≈ val 
end