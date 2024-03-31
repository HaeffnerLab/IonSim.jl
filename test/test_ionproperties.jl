using Test, IonSim
using IonSim.Properties: process_nonlinear_zeeman, loadfromconfig
using Suppressor

@suppress_err begin
    @testset "properties -- Verify nonlinear zeeman shift mapping" begin
        nonlinear_zeeman_config = [
            Dict("level" => "S1/2=0", "sublevel" => 0, "coeffs" => [1.0, 2.0, 3.0]),
            Dict("level" => "S1/2=1", "sublevel" => 0, "coeffs" => [10.0, 20.0, 30.0])
        ]

        processed_functions = process_nonlinear_zeeman(nonlinear_zeeman_config)

        # Really, want to verify that the functions mapped to different 
        # level/sublevel pairs are actually distinct from each other.
        @test isapprox(processed_functions["S1/2=0", 0](2.0), 17.0, rtol=1e-5)
        @test isapprox(processed_functions["S1/2=1", 0](2.0), 170.0, rtol=1e-5)
    end

    @testset "properties -- Verify that all provided ion configs parse" begin
        for ion_fname in readdir("../configs/ions")
            loadfromconfig(joinpath("../configs/ions", ion_fname))
        end
    end
end
