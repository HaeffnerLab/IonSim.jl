using Test, IonSim
using Suppressor
using IonSim:
    linear_chain_normal_modes,
    minimize_psuedopotential_constant,
    minimize_DC_imbalance,
    linear_equilibrium_positions

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
        # test calculation of equilibrium positions for a linear chain of equal mass ions 
        # against known values
        posL = [-2.8708, -2.10003, -1.4504, -0.85378, -0.2821]
        pos = [posL; -1 .* reverse(posL)]
        @test any(isapprox.(IonSim.linear_equilibrium_positions(10), pos, rtol = 1e-4))

        # and test calculation of characterstic length scale for linear chain, equal mass
        C = Ca40()
        @test IonSim.compute_characteristic_length_scale(mass(C), 1e6) ≈ 4.449042804354206e-6

        # and do the same for Anm, which computes the normal modes
        cst = [-0.2132, 0.6742, -0.6742, 0.2132]
        freq, mode = IonSim.linear_chain_normal_modes(
            [1, 1, 1, 1],
            (x = 2, y = 2, z = 1),
            (x = 0, y = 0, z = 1)
        )[2][end]
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

    @testset "ion_configurations -- mixed species" begin
        ####################################################################################
        # mixed-species linear chain tests (axial)
        ####################################################################################

        Ba = 137.9052472 * 1.66054e-27  # mass in kg
        Yb = 170.9363315 * 1.66054e-27

        #= 
        Exact formulas from here:
        Eur. Phys. J. D 13, 261–269 (2001)
        for the axial eigenvalues and eigenvectors of a two-ion mixed-species crystal
        =#
        function two_ion_eigenvalues(k, M)
            m = M[1]
            M = M[2]
            μ = M / m
            C1 = 1 + (1 / μ)
            C2 = sqrt(1 + 1 / μ^2 - 1 / μ)
            return sqrt.((k / m) * [C1 - C2, C1 + C2]) / 2π
        end

        function two_ion_eigenvectors(k, M)
            m = M[1]
            M = M[2]
            μ = M / m
            C = √(1 + μ^2 - μ)
            return [normalize([1 - μ + C, 1] / √μ), normalize([1 - μ - C, 1] / √μ)]
        end

        #=
        Confirm that k_axial = mᵢ(ω_z)ᵢ² ∀i is computed correctly by comparing with 
        value for homogoneous chain
        =#
        M = [Ba]
        com = (x = 3e6, y = 3e6, z = 1e5)

        k_axial_analytic = (M[1] * (2π * com.z)^2)
        k_axial_numeric, a = linear_chain_normal_modes(M, com, ẑ)

        @test k_axial_analytic ≈ k_axial_numeric rtol = 1e-8

        #=
        Compare the axial eigenvalues and eigenvectors for a two-ion mixed-species crystal
        with the analytic formulas.
        =#
        M = [Yb, Ba]
        k_axial, a = linear_chain_normal_modes(M, com, ẑ)

        b, c = two_ion_eigenvalues(k_axial, M)

        @test a[1][1] / a[2][1] ≈ b / c rtol = 1e-8
        @test a[1][1] ≈ b rtol = 1e-8
        @test a[2][1] ≈ c rtol = 1e-8

        vec1, vec2 = two_ion_eigenvectors(k_axial, M)

        # note there is an ambiguity in the overall sign of the vectors
        @test sign(vec1[1]) * vec1 ≈ sign(a[1][2][1]) * a[1][2]
        @test sign(vec2[1]) * vec2 ≈ sign(a[2][2][1]) * a[2][2]

        ####################################################################################
        # mixed-species linear chain tests (radial)
        ####################################################################################

        #=
        Minimization tests. When the linear chain is mixed-species, there is
        no longer generally a unique center-of-mass normal mode. We still allow
        the user to specify a tuple of com values, however. In this case,
        this tuple specifies the lowest (highest) eigenfrequency for the axial (radial)
        direction. In order to find the 
        =#
        M = [Yb, Yb, Yb, Yb, Ba]
        com_x = 4e6 + 1e5 * randn()
        com_x = com_x > 2.8e6 ? com_x : 2.8e6
        com = (x = com_x, y = 3.2e6, z = 1.58e5)

        #= 
        psuedopotential_constant: 
        The psuedopotential scales with mass as α/mᵢ. We call α the psuedopotential_constant. 
        When the user sets com = (x=ωx, y=ωy, z=ωz)/2π and the ions are all the same species, 
        the values com.i correspond to the frequencies
        of the unique com modes. When the ions are different species (so different masses)
        there is no longer generally a unique com mode. In this case we take com to correspond
        to the lowest (highest) eigenfrequency of the axial (radial) modes, which corresponds
        to the homogoneous case. The value com.z is set by the axial spring constant k_axial,
        which is mass-independent. We use com.x to set the psuedopotential constant. Here we 
        check that this is working.
        =#
        pc, a = minimize_psuedopotential_constant(M, com, k_axial)
        @test a[1][1] ≈ com.x rtol = 1e-10

        #= 
        DC_imbalance:
        We perfom the same procedure to set the DC_imbalance based on com.y
        =#
        com_y = com_x + 1e4randn()
        com = (x = com_x, y = com_y, z = 1.58e5)
        res, a = minimize_DC_imbalance(M, com, k_axial, pc)
        @test a[1][1] ≈ com.y rtol = 1e-10

        #=
        Fiducial references: (out of stability regime checks)
        Here we set a fiducial reference based on a normal mode structure that was computed
        at one point (and seemed qualitatively correct) but wasn't strongly cross-checked.
        Note: if one compares our values to this ref: 
        https://doi.org/10.1103/PhysRevA.103.012610
        the structure is qualitatively the same, but there is a difference in the gap between
        the largest eigenfrequency and the rest. Tried to track down the issue, but we seem 
        to be doing evertying correctly on our side.
        =#
        fiducial = [
            (
                3.403e6,
                [
                    5.396022467497479e-5,
                    0.00013583535621031288,
                    0.00043382488271360937,
                    0.002865729612600849,
                    0.9999957890045382
                ]
            )
            (
                2.7468670091375294e6,
                [
                    -0.6308264182921963,
                    -0.5474351959987166,
                    -0.4512008853923827,
                    -0.31430600819209487,
                    0.0014934530345156107
                ]
            )
            (
                2.744229265594708e6,
                [
                    0.6794488606327558,
                    -0.0731278554989956,
                    -0.47502035356991,
                    -0.5543937437967771,
                    0.002191590353902189
                ]
            )
            (
                2.740824881097956e6,
                [
                    -0.35803896744787195,
                    0.6745351565726543,
                    0.12358908202222524,
                    -0.6336653462250121,
                    0.0020947875323122324
                ]
            )
            (
                2.7370687177239084e6,
                [
                    0.11051994247579877,
                    -0.489866035033838,
                    0.7453176887574986,
                    -0.43853920565433885,
                    0.0012320570662960287
                ]
            )
        ]
        com = (x = 3.403e6, y = 3.403e6, z = 1.58e5)
        res, a = minimize_psuedopotential_constant(M, com, k_axial)

        # note there is an ambiguity in the overall sign of the vectors
        for (i, val) in enumerate(a)
            @test sign(val[1][1]) * val[1] ≈ sign(fiducial[i][1][1]) * fiducial[i][1]
            @test sign(val[2][1]) * val[2] ≈ sign(fiducial[i][2][1]) * fiducial[i][2] rtol =
                1e-4
        end

        #=
        For the mixed species case, the user can instead specify the trap parameters in order 
        to calculate the normal mode structure. Here we test that this works and if it 
        corresponds to the com method.
        =#

        #=
        Finally, the user can also specify the individual frequencies and masses of each ion 
        in the chain. We confirm the correspondence of this method with the previous ones.
        =#
    end
end  # end suppress
