using Test, IonSim
using IonSim.PhysicalConstants: m_ca40


@testset "ion_configurations-required_fields" begin
    # linearchain
    C = ca40()
    lc = linearchain(
            ions=[C], com_frequencies=(x=2e6, y=1e6, z=1e6), selected_modes=(x=[],y=[],z=[1])
        )
    @test typeof(label(lc)) == String
    @test typeof(ions(lc)) == Vector{Ion}
end

@testset "ion_configurations-linearchain" begin
    posL = [-2.8708, -2.10003, -1.4504, -0.85378, -0.2821]
    pos = [posL; -1 .* reverse(posL)]
    @test any(isapprox.(linear_equilibrium_positions(10), pos, rtol=1e-4))

    @test characteristic_length_scale(m_ca40, 1e6) ≈ 4.449042804354e-6

    @test_throws AssertionError Anm(2, (x=0.5, y=0.5, z=1), axis=(x=1, y=0, z=0))
    cst = [-0.2132, 0.6742, -0.6742, 0.2132]
    freq, mode = Anm(4, (x=2, y=2, z=1))[end]
    @test freq ≈ √9.308 rtol=1e-4
    if mode[1] > 0
        cst = -1 .* cst
    end
    @test any(isapprox.(mode, cst, rtol=1e-4))

    x = [1e-6, 1e-5]
    IonSim._sparsify!(x, 2e-6)
    @test any(x .== [0, 1e-5])

    C = ca40()
    lc = linearchain(
            ions=[C, C, C, C], com_frequencies=(x=5, y=5, z=1), 
            selected_modes=(x=[], y=[1], z=[4])
        )
    lc.vibrational_modes.y[1].label = "radial"
    lc.vibrational_modes.z[1].label = "axial"
    @test any(isapprox.(lc.vibrational_modes.z[1].mode_structure, mode, rtol=1e-4))
    @test lc["radial"] == lc.vibrational_modes.y[1]

    @test get_vibrational_modes(lc) == [lc["radial"]; lc["axial"]]
end