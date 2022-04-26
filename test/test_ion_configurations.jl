using Test, IonSim
# using IonSim.PhysicalConstants
using Suppressor

@suppress_err begin
    @testset "ion_configurations -- LinearChain" begin
        C = Ca40()
        lc = LinearChain(
            ions = [C, C, C, C],
            com_frequencies = (x = 5, y = 5, z = 1),
            vibrational_modes = (y = [1], z = [4])
        )
        @test ions(lc) == lc.ions
        # test get_vibrational_modes, which should return an array of the selected
        # VibrationalModes in the linear chain
        vms = lc.vibrational_modes
        @test get_vibrational_modes(lc) == [vms.x..., vms.y..., vms.z...]

        # make sure ion numbers are updated
        @test [ionnumber(I) for I in lc.ions] == [1, 2, 3, 4]

        # cmode -> pre-evaluated equilibrium positions for four ion chain
        cmode = [-0.06392393573, -0.0202155287427, 0.0202155287427, 0.0639239357]
        if ionposition(lc.ions[1]) > 0
            cmode .*= -1
        end
        @test [ionposition(I) for I in lc.ions] ≈ cmode rtol = 1e-6

        # should get warning if same ion is input multiple times to ion kwarg
        warning = "Some ions point to the same thing. Making copies."
        @test_logs (:warn, warning) LinearChain(
            ions = [C, C],
            com_frequencies = (x = 4, y = 4, z = 1),
            vibrational_modes = (x = [], y = [], z = [1, 2])
        )
        # and copies should be made of the repeated ions, so that they are no longer the same
        chain1 = LinearChain(
            ions = [C, C],
            com_frequencies = (x = 4, y = 4, z = 1),
            vibrational_modes = (x = [], y = [], z = [1, 2])
        )
        @test !(chain1.ions[1] ≡ chain1.ions[2])

        # Make sure there are no errors with print/show
        @suppress print(lc)
        show(lc)
    end

    @testset "ion_configurations -- general" begin
        # test calculation of equilibrium positions for a linear chain of equal mass ions against
        # known values
        posL = [-2.8708, -2.10003, -1.4504, -0.85378, -0.2821]
        pos = [posL; -1 .* reverse(posL)]
        @test any(isapprox.(linear_equilibrium_positions(10), pos, rtol = 1e-4))

        # and test calculation of characterstic length scale for linear chain, equal mass
        C = Ca40()
        @test characteristic_length_scale(mass(C), 1e6) ≈ 4.449042804354206e-6

        # and do the same for Anm, which computes the normal modes
        @test_throws AssertionError Anm(2, (x = 0.5, y = 0.5, z = 1), (x = 1, y = 0, z = 0))
        cst = [-0.2132, 0.6742, -0.6742, 0.2132]
        freq, mode = Anm(4, (x = 2, y = 2, z = 1), (x = 0, y = 0, z = 1))[end]
        @test freq ≈ √9.308 rtol = 1e-4
        if mode[1] > 0
            cst = -1 .* cst
        end
        @test any(isapprox.(mode, cst, rtol = 1e-4))

        # test '_sparsify', which should drop small values from an array
        x = [1e-6, 1e-5]
        IonSim._sparsify!(x, 2e-6)
        @test any(x .== [0, 1e-5])
    end
end  # end suppress
