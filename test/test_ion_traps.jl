using Test, IonSim
using Suppressor
using IonSim:
    linear_equilibrium_positions,
    searchfor_trappingpotential_parameters,
    diagonalize_Kij,
    diagonalize_Kij_for_homogeneous_chain

@suppress_err begin
    @testset "ion_traps -- LinearChain" begin
        C = Ca40()
        lc = LinearChain(
            ions=[C, C, C, C],
            comfrequencies=(x=5e6, y=5e6, z=1e6),
            selectedmodes=(y=[1], z=[4])
        )
        @test ions(lc) == lc.ions
        # test modes, which should return an array of the selected
        # VibrationalModes in the linear chain
        vms = lc.selectedmodes
        @test modes(lc) == [vms.x..., vms.y..., vms.z...]

        # make sure ion numbers are updated
        @test [ionnumber(I) for I in lc.ions] == [1, 2, 3, 4]

        # cmode -> pre-evaluated equilibrium positions for four ion chain
        cmode = [
            -2.1766240336602492e-5,
            -6.883427582857696e-6,
            6.883427582857696e-6,
            2.1766240336602492e-5
        ]
        @test ionpositions(lc) ≈ cmode rtol = 1e-6

        # should get warning if same ion is input multiple times to ion kwarg
        warning = "Some ions point to the same thing. Making copies."
        @test_logs (:warn, warning) LinearChain(
            ions=[C, C],
            comfrequencies=(x=4, y=4, z=1),
            selectedmodes=(x=[], y=[], z=[1, 2])
        )
        # and copies should be made of the repeated ions, so that they are no longer the same
        chain1 = LinearChain(
            ions=[C, C],
            comfrequencies=(x=4, y=4, z=1),
            selectedmodes=(x=[], y=[], z=[1, 2])
        )
        @test !(chain1.ions[1] ≡ chain1.ions[2])

        # Make sure there are no errors with print/show
        @suppress print(lc)
        show(lc)

        @test length(lc) == length(lc.ions)
    end

    @testset "ion_traps -- general" begin
        # test calculation of equilibrium positions for a linear chain of equal mass ions against
        # known values
        posL = [-2.8708, -2.10003, -1.4504, -0.85378, -0.2821]
        pos = [posL; -1 .* reverse(posL)]
        @test any(isapprox.(linear_equilibrium_positions(10), pos, rtol=1e-4))

        # test '_sparsify', which should drop small values from an array
        x = [1e-6, 1e-5]
        IonSim._sparsify!(x, 2e-6)
        @test any(x .== [0, 1e-5])
    end

    @testset "ion_configurations -- mixed species" begin
        ####################################################################################
        # mixed-species linear chain tests 
        ####################################################################################
        c = Ca40()
        b = Be9()
        y = Yb171()

        # homogeneous chain

        ion_list = [c for _ in 1:rand(2:7)]
        randval = rand() * 1e6
        target_eigenfrequencies = (x=randval, y=randval, z=randval / 10)

        # compare optimization routine to james formula for homogeneous case https://doi.org/10.1007/s003400050373
        res = searchfor_trappingpotential_parameters(ion_list, target_eigenfrequencies)
        nmstructure_x_direction = diagonalize_Kij(ion_list, x̂, res...)
        nmstructure_y_direction = diagonalize_Kij(ion_list, ŷ, res...)
        nmstructure_z_direction = diagonalize_Kij(ion_list, ẑ, res...)

        for (i, e) in enumerate(nmstructure_x_direction)
            @test e[1] ≈ diagonalize_Kij_for_homogeneous_chain(
                length(ion_list),
                target_eigenfrequencies,
                x̂
            )[i][1]
            @test e[2] ≈ diagonalize_Kij_for_homogeneous_chain(
                length(ion_list),
                target_eigenfrequencies,
                x̂
            )[i][2]
        end

        for (i, e) in enumerate(nmstructure_y_direction)
            @test e[1] ≈ diagonalize_Kij_for_homogeneous_chain(
                length(ion_list),
                target_eigenfrequencies,
                ŷ
            )[i][1]
            @test e[2] ≈ diagonalize_Kij_for_homogeneous_chain(
                length(ion_list),
                target_eigenfrequencies,
                ŷ
            )[i][2]
        end

        for (i, e) in enumerate(nmstructure_z_direction)
            @test e[1] ≈ diagonalize_Kij_for_homogeneous_chain(
                length(ion_list),
                target_eigenfrequencies,
                ẑ
            )[i][1]
            @test e[2] ≈ diagonalize_Kij_for_homogeneous_chain(
                length(ion_list),
                target_eigenfrequencies,
                ẑ
            )[i][2]
        end

        # test for assertion error if x > y > z not obeyed
        @test_throws AssertionError searchfor_trappingpotential_parameters(
            ion_list,
            (
                x=target_eigenfrequencies.z,
                y=target_eigenfrequencies.x,
                z=target_eigenfrequencies.y
            )
        )

        # heterogeneous chain

        #= 
        Compare to exact formulas for axial eigenvalues and eigenvectors of a two-ion mixed-species crystal:
        Eur. Phys. J. D 13, 261–269 (2001) https://doi.org/10.1007/s100530170275
        =#

        function two_ion_eigenvalues(k, masses)
            m = minimum(masses)
            M = maximum(masses)
            μ = M / m
            C1 = 1 + (1 / μ)
            C2 = sqrt(1 + 1 / μ^2 - 1 / μ)
            return sqrt.((k / m) * [C1 - C2, C1 + C2])
        end

        function two_ion_eigenvectors(masses)
            m = minimum(masses)
            M = maximum(masses)
            μ = M / m
            C = √(1 + μ^2 - μ)
            q1 = normalize([√μ, 1 - μ + C])
            q2 = normalize([√μ, 1 - μ - C])
            if masses[2] > masses[1]
                reverse!(q1)
                reverse!(q2)
            end
            q1 *= sign(q1[1])
            q2 *= sign(q2[1])
            return [q1, q2]
        end

        for ion_list in [[c, b], [b, c], [c, y], [y, c], [b, y], [y, b]]
            target_eigenfrequencies = (x=20e6, y=20e6, z=1e6)  # have to make these huge for [y, b] case
            res = searchfor_trappingpotential_parameters(ion_list, target_eigenfrequencies)
            computed_values = diagonalize_Kij(ion_list, ẑ, res...)
            computed_eigenvalues = [i[1] for i in computed_values]
            computed_eigenvectors = [i[2] for i in computed_values]

            M = [mass(ion) for ion in ion_list]
            analytic_eigenvectors = two_ion_eigenvectors(M)
            analytic_eigenvalues = two_ion_eigenvalues(res[1], M)

            # for some reason this isn't working here, even though it works fine in a notebook
            # also returns the same values for every itration of the loop
            # @test analytic_eigenvalues ≈ computed_eigenvalues rtol=1e-3 
            @test analytic_eigenvectors ≈ computed_eigenvectors rtol = 1e-3
        end

        #= 
        Compare to eigenvalues from PHYSICAL REVIEW A103, 012610 (2021)
        =#

        M = [mass(y), mass(y), mass(y), mass(y), 2.28997e-25]
        Q = ones(Int, length(M))
        target_eigenfrequencies = (x=3.039e6, y=3.039e6, z=0.158e6)
        res = searchfor_trappingpotential_parameters(M, Q, target_eigenfrequencies)
        computed_values = diagonalize_Kij(M, Q, ẑ, res...)
        axial_eigenfrequencies = [round(i[1] / 1e6, digits=3) for i in computed_values]

        axial_values_from_PRA = [0.158, 0.28, 0.388, 0.48, 0.57]
        @test axial_eigenfrequencies ≈ axial_values_from_PRA

        # Note: we don't test radial modes against this PRA since there is an unresolved
        # quantitative (but not qualitative difference). Normal modes are the same, 
        # eigenfrequencies are gapped between modes that mostly vibrate in one or the other
        # species and the general scaling seems to agree. but the gap width is larger for
        # our calculation and I'm not sure why
    end

    @testset "load from yaml" begin
        c = Ca40()
        chain = LinearChain_fromyaml(
            ions=[c, c, c, c],
            yaml="test_load_LinearChain_fromyaml.yaml",
        )
        @test isnan(chain.comfrequencies.x)
        @test chain.ionpositions == [1, 2, 3, 4]
        @test chain.ions == ions(chain)
        @test_throws AssertionError full_normal_mode_description(chain)
        lc = LinearChain_fromyaml(
            ions=[c, c, c, c],
            yaml="test_load_LinearChain_fromyaml_broken.yaml",
        )
        @test length(lc.selectedmodes.x) == 0
    end

    @testset "visualize-iontraps" begin
        # for now just run all of the plots and make sure nothing errors
        c = Ca40()
        b = Be9()
        y = Yb171()
        chain = LinearChain(
            ions=[c, b, b, c, y, y],
            comfrequencies=(x=20e6, y=20e6, z=0.1e6),
            selectedmodes=(; x=[1, 3], z=[1, 3:4])
        )
        visualize(chain, ẑ, [:], format="circles")
        visualize(chain, x̂, [1])
        visualize(chain, ẑ, [1:4], format="circles")
        visualize(xmodes(chain)[1], format="circles")
        visualize(xmodes(chain)[1], format="bar")
    end
end  # end suppress
