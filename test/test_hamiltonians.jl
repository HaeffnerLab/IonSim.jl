using QuantumOptics: DenseOperator
using LinearAlgebra: diagm, norm
using Test, IonSim
using Suppressor
using Unitful
using IonSim.PhysicalConstants: INVERSE_TIME #todo: i shouldn't have to do this, it should be exported

Base.:(*)(s::String, i::Int) = i == 1 ? s : 0
Base.:(*)(i::Int, s::String) = Base.:(*)(s, i)
Base.:(+)(s::String, i::Int) = i == 0 ? s : nothing
Base.:(+)(i::Int, s::String) = Base.:(+)(s, i)

function get_indices(
    ion_participation,
    vibrational_mode_dimension,
    two_modes;
    cutoff = [Inf, Inf]
)
    N = vibrational_mode_dimension[1]
    σ1 = [0 0; "a" 0]
    σ2 = [0 0; "b" 0]
    I = [1 0; 0 1]
    c = ["c$i$j," for i in 1:N, j in 1:N]
    if two_modes
        M = vibrational_mode_dimension[2]
        d = ["d$i$j," for i in 1:M, j in 1:M]
        subarrays = [[d, c, I, σ1], [d, c, σ2, I]]
    else
        subarrays = [[c, I, σ1], [c, σ2, I]]
    end
    if ion_participation == 1
        k = kron(subarrays[1]...)
    elseif ion_participation == 2
        k = kron(subarrays[2]...)
    elseif ion_participation == 12
        k = kron(subarrays[1]...) + kron(subarrays[2]...)
    end
    ridxs = Dict()
    cidxs = Dict()
    nzidxs = []
    for i in 1:size(k, 2), j in 1:size(k, 2)
        if k[i, j] == 0
            continue
        end
        mn =
            map(x -> [parse(Int, x[2]), parse(Int, x[3])], split(k[i, j], ",")[1:(end - 1)])
        if sum([abs(x[1] - x[2]) > cutoff[l] for (l, x) in enumerate(mn)]) > 0
            continue
        end
        if i > j
            if !(k[i, j] in keys(ridxs))
                ridxs[k[i, j]] = [(i, j)]
            else
                push!(ridxs[k[i, j]], (i, j))
            end
        else
            if !(k[i, j] in keys(cidxs))
                parity = sum(map(x -> isodd(abs(x[1] - x[2])), mn))
                cidxs[k[i, j]] = [(-1 * isodd(parity), 0), (i, j)]
            else
                push!(cidxs[k[i, j]], (i, j))
            end
        end
    end
    return collect(values(ridxs)), collect(values(cidxs))
end

@suppress_err begin
    @testset "hamiltonians -- misc functions" begin
        # TODO: enforce units on here
        # make sure _D and _D_constant_eta return the same thing for constant eta
        L = 10
        Ω = 1e6randn()u"1/s"
        Δ = 10randn()u"1/s"
        η = [abs(randn()) for _ in 1:L]
        ν = [10randn() for _ in 1:L]u"1/s"
        timescale = 10randn()u"s"
        n = [[rand(1:20) for _ in 1:L], [rand(1:5) for _ in 1:L]]
        t = randn()
        d1 = IonSim._D(Ω, Δ, η, ν, timescale, n, t, L)
        D = [IonSim._Dnm_cnst_eta(η[i], n[1][i], n[2][i]) for i in 1:L]
        d2 = IonSim._D_cnst_eta(Ω, Δ, ν, timescale, n, D, t, L)
        @test abs(sum(d1 .- d2)) < 1e-8

        # _flattenall (this function should completely flatten nested arrays)
        a = [[[1], [2], [3]], [[4], [5], [6]]]
        @test IonSim._flattenall(a) == collect(1:6)

        # _get_kron_indxs
        a = ["a$i$j," for i in 1:3, j in 1:3]
        b = ["b$i$j," for i in 1:3, j in 1:3]
        c = ["c$i$j," for i in 1:3, j in 1:3]
        abc = kron(a, b, c)
        z = 1:3
        i1 = rand(z)
        i2 = rand(z)
        j1 = rand(z)
        j2 = rand(z)
        k1 = rand(z)
        k2 = rand(z)
        ijk = [(i1, i2), (j1, j2), (k1, k2)]
        (I, J) = IonSim._get_kron_indxs([(i1, i2), (j1, j2), (k1, k2)], [3, 3, 3])
        @test [
            (parse(Int, i[2]), parse(Int, i[3])) for i in split(abc[I, J], ",")[1:(end - 1)]
        ] == ijk
        # ^this tests that  abc[I, J] = a[i1,i2] * a[j1,j2] * c[k1,k2]

        # _inv_get_kron_indxs
        indxs = IonSim._inv_get_kron_indxs((I, J), (3, 3, 3))
        indxs = collect(zip(indxs[1], indxs[2]))
        @test indxs == ijk
    end

    @testset "hamiltonians -- parameter arrays" begin
        # setup system
        C = Ca40([("S1/2", -1 / 2), ("D5/2", -1 / 2)])
        C1 = copy(C)

        L1 = Laser()
        L1.pointing = [(1, 1.0), (2, 1.0)]
        L1.λ = transitionwavelength(C, ("S1/2", "D5/2"))
        L2 = copy(L1)

        # _Ωmatrix
        L1.E = 1u"V/m"
        L1.ϕ = 0
        L2.E = 2u"V/m"
        L2.ϕ = 2
        chain = LinearChain(
            ions = [C, C1],
            com_frequencies = (x = 3e6u"1/s", y = 3e6u"1/s", z = 1e6u"1/s"),
            vibrational_modes = (; z = [1])
        )
        T = Trap(configuration = chain, lasers = [L1, L2])
        Ωnmkj = IonSim._Ωmatrix(T, 1u"s")
        t = rand(0:1e-3:100)

        # coupling strength between ion1-laser1 and ion2-laser2 should be identical for all
        # transitions since the ions are identical
        resolve(x, t) = typeof(x) <: Number ? x : x(t)
        @test [resolve(i, t) for i in Ωnmkj[1, 1]] == [resolve(i, t) for i in Ωnmkj[2, 1]]

        # coupling strength between ion1-laser1 and ion1-laser2 should be proportional according
        # to the factor 2exp(-2im) due to the differences we've set in L2.E and L2.ϕ
        @test [2exp(-2im) * resolve(i, t) for i in Ωnmkj[1, 1]] == [resolve(i, t) for i in Ωnmkj[1, 2]]
        @test [2exp(-2im) * resolve(i, t) for i in Ωnmkj[2, 1]] == [resolve(i, t) for i in Ωnmkj[2, 2]]

        # make sure time-dep L.E and L.ϕ propagate appropriately
        L1.E = x -> 1u"V/m" * cos(x)
        L1.ϕ = t -> t^2
        Ωnmkj = IonSim._Ωmatrix(T, 1u"s")
        t = 0:1e-3:100
        for Ω in Ωnmkj[1, 1]
            if typeof(Ω) <: Number
                continue
            end
            @test Ω.(t) ≈ @.(Ω(0) * cos(t) * exp(-im * t^2))
        end

        # _Δmatrix
        # zero B, zero laser detuning, zero stark shift should give an array of zeros
        Δ = IonSim._Δmatrix(T, 1u"s")
        for i in 1:2, j in 1:2
            @test Δ[i, j][1] == 0
        end

        # zero B, zero stark shift, non-zero laser detuning should return columns that all have
        # the same value
        L1.Δ = 1u"1/s"
        L2.Δ = -1u"1/s"
        Δ = IonSim._Δmatrix(T, 1u"s")
        @test Δ[1, 1][1] ≈ 2π && Δ[2, 1][1] ≈ 2π && Δ[1, 2][1] ≈ -2π && Δ[2, 2][1] ≈ -2π

        # zero B, zero laser detuning, now add stark shift to just one of the ions
        L1.Δ = 0u"1/s"
        L2.Δ = 0u"1/s"
        set_stark_shift!(C, ("S1/2", -1 / 2), 1u"1/s")
        Δ = IonSim._Δmatrix(T, 1u"s")
        @test Δ[1, 1][1] ≈ 2π && Δ[1, 2][1] ≈ 2π && Δ[2, 1][1] ≈ 0 && Δ[2, 2][1] ≈ 0

        # lastly let's test when resonant
        zero_stark_shift!(C)
        T.B = 1e-4u"T"
        L1.λ = transitionwavelength(C, (("S1/2", -1 / 2), ("D5/2", -1 / 2)), T)
        Δ = IonSim._Δmatrix(T, 1u"s")
        @test Δ[1, 1][1] ≈ 0

        # _ηmatrix
        chain = LinearChain(
            ions = [C, C1],
            com_frequencies = (x = 2e6u"1/s", y = 2e6u"1/s", z = 1e6u"1/s"),
            vibrational_modes = (x = [1], y = [2], z = [1])
        )
        L1.k = (x̂ + ẑ) / √2
        L2.k = (ŷ + ẑ) / √2
        L3 = Laser()
        L3.pointing = [(1, 1.0), (2, 1.0)]
        L3.λ = transitionwavelength(C, ("S1/2", "D5/2"))
        L3.k = ẑ
        T = Trap(configuration = chain, lasers = [L1, L2, L3])
        η = IonSim._ηmatrix(T)

        # η[1,1,1] corresponds laser1, ion1, mode1 (x̂) and η[1,1,3] corresponds laser1, ion1,
        # mode3 (ẑ). They both have the same projection on L1.k, but the frequency of mode1 is
        # twice as large as mode2, so it's L-D factor should be √2 times smaller
        @test abs(η[1, 1, end](1.0)) ≈ abs(η[1, 1, end - 2](1.0) / √2)
        # With this setup, the L-D factor should be the same for ion1 and ion2
        @test η[1, 1, end](1.0) ≈ η[2, 1, end](1.0)
        # η[1, 2, 2] is the 1st ion, 2nd laser and y-stretch-mode. The L-D factor should be
        # opposite in sign, equal in magnitude to the 2nd ion, 2nd laser and y-stretch-mode
        @test η[1, 2, end - 1](1.0) ≈ -η[2, 2, end - 1](1.0)
        # L3, which is in the ẑ direction should only have projection on zmode (mode3)
        @test η[1, 3, end] ≡ 0
        @test η[1, 3, end - 1] ≡ 0
        @test η[1, 3, end - 2](0.0) != 0
        # test construction of time-dep δν. If δν = 1e6*t, then after 3e-6 seconds (and since
        # ν=1e6), √(ν+δν(t)) = √2 * √(ν) = √2 * √(ν+δν(0))
        chain.vibrational_modes.z[1].δν = t -> 1e6t
        η = IonSim._ηmatrix(T)
        @test η[1, 1, end - 2](0.0) ≈ 2 * η[1, 1, end - 2](3)
    end

    @testset "hamiltonians -- setup functions" begin
        # tests for _setup_global_B_hamiltonian

        # setup system
        C = Ca40()
        L = Laser()
        L.λ = transitionwavelength(C, ("S1/2", "D5/2"))
        chain = LinearChain(
            ions = [C],
            com_frequencies = (x = 3e6u"1/s", y = 3e6u"1/s", z = 1e6u"1/s"),
            vibrational_modes = (; z = [1])
        )
        T = Trap(configuration = chain, lasers = [L], δB = 0u"T")
        global_B_indices, global_B_scales, bfunc =
            IonSim._setup_global_B_hamiltonian(T, 1u"s")

        # T.δB = 0 -> global_B_indices and global_B_scales should be empty arrays
        @test length(global_B_indices) == 0 && length(global_B_scales) == 0
        # and bfunc should be identically zero like T.δB
        @test sum(bfunc.(0:1e4)) == 0

        # now let's test nontrivial T.δB
        T.δB = x -> 1u"T" * sin(x)
        t = 0:0.1:10
        T.configuration.vibrational_modes.z[1].N = 3
        global_B_indices, global_B_scales, bfunc =
            IonSim._setup_global_B_hamiltonian(T, 1u"s")
        # make sure bfunc is correct
        @test bfunc.(t) == 2π .* sin.(t)
        # make sure all energy levels are being recorded
        @test sort(IonSim._flattenall(global_B_indices)) == collect(1:32)
        # Not sure what this test is supposed to be; commented out for now
        # # and make sure that the susceptibilites are correct
        # zs = [zeeman_shift(1, quantumnumbers(C, sublevel)) for sublevel in sublevels(C)]
        # @test length(unique([global_B_scales; zs])) == length(zs)

        # test _setup_δν_hamiltonian

        # setup system
        chain = LinearChain(
            ions = [C],
            com_frequencies = (x = 3e6u"1/s", y = 3e6u"1/s", z = 1e6u"1/s"),
            vibrational_modes = (; z = [1])
        )
        T = Trap(configuration = chain, lasers = [L], δB = 0u"T")
        # should return empty arrays if δν=0
        δν_indices, δν_functions = IonSim._setup_δν_hamiltonian(T, 1u"s")
        @test length(δν_indices) == 0 && length(δν_functions) == 0

        # test output for simple case
        T.configuration.vibrational_modes.z[1].N = 3
        T.configuration.vibrational_modes.z[1].δν = 1
        δν_indices, δν_functions = IonSim._setup_δν_hamiltonian(T, 1u"s")
        indxs = [[8j + i for i in 1:8] for j in 1:3]
        @test δν_indices[1] == indxs

        # add another mode with δν=0 and test output
        chain = LinearChain(
            ions = [C],
            com_frequencies = (x = 3e6u"1/s", y = 3e6u"1/s", z = 1e6u"1/s"),
            vibrational_modes = (y = [1], z = [1])
        )
        T = Trap(configuration = chain, lasers = [L], δB = 0u"T")
        T.configuration.vibrational_modes.y[1].N = 3
        T.configuration.vibrational_modes.z[1].N = 3
        T.configuration.vibrational_modes.z[1].δν = 1
        δν_indices, δν_functions = IonSim._setup_δν_hamiltonian(T, 1u"s")
        indxs1 = [[8 * 4 * j + i for i in 1:(8 * 4)] for j in 1:3]
        @test δν_indices[1] == indxs1

        # test output when both modes have nonzero δν
        T.configuration.vibrational_modes.z[1].δν = 1
        T.configuration.vibrational_modes.y[1].δν = 1
        δν_indices, δν_functions = IonSim._setup_δν_hamiltonian(T, 1u"s")
        indxs = [[8 * j + i for i in 1:(8 * 8)] for j in 1:3]
        @test δν_indices[1][1] ==
              [vcat([[8 * (4(j - 1) + 1) + i for i in 1:8] for j in 1:4]...);]
        @test δν_indices[2] == indxs1
        # we have two modes but only one nontrivial δν, so length(δν_functions) should equal 1
        @test length(δν_functions) == 2

        # finally, make sure δν_functions are being constructed appropriately
        T.configuration.vibrational_modes.y[1].δν = cos
        T.configuration.vibrational_modes.z[1].δν = sin
        δν_indices, δν_functions = IonSim._setup_δν_hamiltonian(T, 1u"s")
        t = 0:0.1:10
        @test δν_functions[1].(t) == 2π .* cos.(t)
        @test δν_functions[2].(t) == 2π .* sin.(t)

        # _setup_fluctuation_hamiltonian
        T = Trap(configuration = chain, lasers = [L], δB = 1u"T")
        T.configuration.vibrational_modes.y[1].δν = cos
        T.configuration.vibrational_modes.z[1].δν = sin
        δν_indices, δν_functions = IonSim._setup_δν_hamiltonian(T, 1u"s")
        global_B_indices, global_B_scales, bfunc =
            IonSim._setup_global_B_hamiltonian(T, 1u"s")
        all_unique_indices, gbi, gbs, bfunc1, δνi, δνfuncs =
            IonSim._setup_fluctuation_hamiltonian(T, 1u"s")
        # make sure information from _setup_global_B_hamiltonian and _setup_δν_hamiltonian
        # is propagated appropriately
        @test any(gbi .== global_B_indices)
        @test any(gbs .== global_B_scales)
        @test bfunc(17.0) == bfunc1(17.0)
        @test any(δνi .== δν_indices)
        @test δνfuncs[end](92.0) == δν_functions[end](92.0)

        # _setup_base_hamiltonian
        # Let's make sure that _setup_base_hamiltonian is recording the appropriate indices.
        C = Ca40([("S1/2", -1 / 2), ("D5/2", -1 / 2)])
        L = Laser()
        L.λ = transitionwavelength(C, ("S1/2", "D5/2"))
        chain = LinearChain(
            ions = [C, C],
            com_frequencies = (x = 3e6u"1/s", y = 3e6u"1/s", z = 1e6u"1/s"),
            vibrational_modes = (; z = [1])
        )
        T = Trap(
            configuration = chain,
            B = 4e-4u"T",
            Bhat = (x̂ + ŷ + ẑ) / √3,
            lasers = [L]
        )
        mode = T.configuration.vibrational_modes.z[1]
        mode.N = rand(1:8)
        N = mode.N + 1
        Efield_from_rabi_frequency!(
            1e6u"1/s",
            T,
            1,
            1,
            (("S1/2", -1 / 2), ("D5/2", -1 / 2))
        )

        ## first just shine light on 1st ion
        L.pointing = [(1, 1.0)]
        ridxs, cidxs = get_indices(1, [N], false)
        _, r, c = IonSim._setup_base_hamiltonian(T, 1e-6u"s", 100, Inf, "analytic", true)
        @test length(unique([ridxs; r])) == (N^2 + N) / 2 == length(ridxs)
        @test length(unique([cidxs; c])) - 1 == (N^2 - N) / 2 == length(cidxs)

        ## now just shine light on second ion
        L.pointing = [(2, 1.0)]
        ridxs, cidxs = get_indices(2, [N], false)
        _, r, c = IonSim._setup_base_hamiltonian(T, 1e-6u"s", 100, Inf, "analytic", true)
        @test length(unique([ridxs; r])) == (N^2 + N) / 2 == length(ridxs)
        @test length(unique([cidxs; c])) - 1 == (N^2 - N) / 2 == length(cidxs)

        ## now shine light on both ions
        L.pointing = [(1, 1.0), (2, 1.0)]
        ridxs, cidxs = get_indices(12, [N], false)
        _, r, c = IonSim._setup_base_hamiltonian(T, 1e-6u"s", 100, Inf, "analytic", true)
        @test length(unique([ridxs; r])) == (N^2 + N) == length(ridxs)
        @test length(unique([cidxs; c])) - 1 == (N^2 - N) == length(cidxs)

        ## test for two modes
        C1 = Ca40([("S1/2", -1 / 2), ("D5/2", -1 / 2)])
        L1 = Laser()
        L1.λ = transitionwavelength(C1, ("S1/2", "D5/2"))
        chain1 = LinearChain(
            ions = [C1, C1],
            com_frequencies = (x = 3e6u"1/s", y = 3e6u"1/s", z = 1e6u"1/s"),
            vibrational_modes = (; z = [1, 2])
        )
        T1 = Trap(
            configuration = chain1,
            B = 4e-4u"T",
            Bhat = (x̂ + ŷ + ẑ) / √3,
            lasers = [L1]
        )
        mode1 = T1.configuration.vibrational_modes.z[1]
        mode2 = T1.configuration.vibrational_modes.z[2]
        mode1.N = N - 1
        mode2.N = rand(1:8)
        M = mode2.N + 1
        NM = N * M
        Efield_from_rabi_frequency!(
            1e6u"1/s",
            T1,
            1,
            1,
            (("S1/2", -1 / 2), ("D5/2", -1 / 2))
        )

        ridxs, cidxs = get_indices(12, [N, M], true)
        _, r, c = IonSim._setup_base_hamiltonian(T, 1e-6u"s", 100, Inf, "analytic", true)
        @test length(unique([ridxs; r])) == (NM^2 + NM) == length(ridxs)
        @test length(unique([cidxs; c])) - 1 == (NM^2 - NM) == length(cidxs)

        ## now set non-trivial lamb_dicke_order
        lamb_dicke_order = [rand(1:N), rand(1:M)]
        ridxs, cidxs = get_indices(12, [N, M], true, cutoff = reverse(lamb_dicke_order))
        _, r, c = IonSim._setup_base_hamiltonian(
            T1,
            1e-6u"s",
            lamb_dicke_order,
            Inf,
            "analytic",
            true
        )
        @test length(unique([ridxs; r])) == length(ridxs)
        @test length(unique([cidxs; c])) - 1 == length(cidxs)

        ## when laser tuned to carrier, setting an rwa_cutoff below the vibrationl_mode frequency
        ## should have the same effect
        L.λ = transitionwavelength(C1, (("S1/2", -1 / 2), ("D5/2", -1 / 2)), T1)
        _, repeated_indices, conj_repeated_indices =
            IonSim._setup_base_hamiltonian(T1, 1u"s", 100, 1e5, "analytic", true)
        @test length(unique([ridxs; r])) == length(ridxs)
        @test length(unique([cidxs; c])) - 1 == length(cidxs)
    end

    # @testset "hamiltonian" begin
    #     C1 = Ca40([("S1/2", -1/2), ("D5/2", -1/2)])
    #     C2 = Ca40([("S1/2", -1/2), ("D5/2", -1/2)])
    #     L = Laser()
    #     L.pointing = [(1, 1.0), (2, 1.0)]
    #     L.λ = transitionwavelength(C1, ("S1/2", "D5/2"))
    #     chain = LinearChain(
    #                 ions=[C1, C2], com_frequencies=(x=3e6,y=3e6,z=1e6),
    #                 vibrational_modes=(;z=[1,2])
    #             )
    #     T = Trap(configuration=chain, B=4e-4, Bhat=(x̂ + ŷ + ẑ)/√3, lasers=[L])
    #     mode1 = T.configuration.vibrational_modes.z[1]
    #     mode2 = T.configuration.vibrational_modes.z[2]
    #     Δ = round(randn(), digits=5) * 1e5  # TODO: this begins to fail at below 1 Hz!
    #     L.Δ =  Δ
    #     ϕ = randn()
    #     L.ϕ = ϕ
    #     mode1.N = 10
    #     mode2.N = 9
    #     Ω = randn()
    #     Efield_from_rabi_frequency!(Ω * 1e6, T, 1, 1, (("S1/2", -1/2), ("D5/2", -1/2)))

    #     # full hamiltonian (w conj_repeated_indices)
    #     timescale = 1e-6
    #     ion_op(t) = Ω * π * exp(-im * (2π * Δ*t*timescale + ϕ)) * C1[("D5/2", -1/2)] ⊗ C1[("S1/2", -1/2)]'
    #     η11 = get_η(mode1, L, C1)
    #     η12 = get_η(mode2, L, C1)
    #     η21 = get_η(mode1, L, C2)
    #     η22 = get_η(mode2, L, C2)
    #     mode_op1(t; η) = displace(mode1, im * η * exp(im * 2π * t), method="truncated")
    #     mode_op2(t; η) = displace(mode2, im * η * exp(im * 2π * √3 * t), method="truncated")
    #     Hp(t) = (ion_op(t) ⊗ one(C2) ⊗ mode_op1(t, η=η11) ⊗ mode_op2(t, η=η12)
    #              + one(C1) ⊗ ion_op(t) ⊗ mode_op1(t, η=η21) ⊗ mode_op2(t, η=η22))
    #     qoH(t) = Hp(t) + dagger(Hp(t))
    #     tp = abs(51randn())

    #     H = hamiltonian(T, lamb_dicke_order=101)
    #     H1 = hamiltonian(T, lamb_dicke_order=101, time_dependent_eta=true)
    #     @test norm(qoH(tp).data - H(tp, 0).data) < 1e-4
    #     @test norm(qoH(tp).data - H1(tp, 0).data) < 1e-4

    #     mode_op1(t; η) = displace(mode1, im * η * exp(im * 2π * t), method="analytic")
    #     mode_op2(t; η) = displace(mode2, im * η * exp(im * 2π * √3 * t), method="analytic")
    #     H1 = hamiltonian(T, lamb_dicke_order=101, displacement="analytic", time_dependent_eta=false)
    #     H = hamiltonian(T, lamb_dicke_order=101, displacement="analytic", time_dependent_eta=true)
    #     @test norm(qoH(tp).data - H(tp, 0).data) < 1e-4
    #     @test H1(tp, 0).data ≈ H(tp, 0).data

    #     # full hamiltonian (w/o conj_repeated_indices)
    #     H = hamiltonian(T, lamb_dicke_order=101, rwa_cutoff=1e10)
    #     H1 = hamiltonian(T, lamb_dicke_order=101, time_dependent_eta=true, rwa_cutoff=1e10)
    #     mode_op1(t; η) = displace(mode1, im * η * exp(im * 2π * t), method="truncated")
    #     mode_op2(t; η) = displace(mode2, im * η * exp(im * 2π * √3 * t), method="truncated")
    #     @test norm(qoH(tp).data - H(tp, 0).data) < 1e-4
    #     @test norm(qoH(tp).data - H1(tp, 0).data) < 1e-4

    #     mode_op1(t; η) = displace(mode1, im * η * exp(im * 2π * t), method="analytic")
    #     mode_op2(t; η) = displace(mode2, im * η * exp(im * 2π * √3 * t), method="analytic")
    #     H1 = hamiltonian(T, lamb_dicke_order=101, displacement="analytic", time_dependent_eta=false, rwa_cutoff=1e10)
    #     H = hamiltonian(T, lamb_dicke_order=101, displacement="analytic", time_dependent_eta=true, rwa_cutoff=1e10)
    #     @test norm(qoH(tp).data - H(tp, 0).data) < 1e-4
    #     @test H1(tp, 0).data ≈ H(tp, 0).data

    #     # # Lamb-Dicke
    #     mode_op11 = DenseOperator(mode1, diagm(0 => [1 - η11^2 * i for i in 0:mode1.N]))
    #     mode_op12 = DenseOperator(mode2, diagm(0 => [1 - η12^2 * i for i in 0:mode2.N]))
    #     mode_op21 = DenseOperator(mode1, diagm(0 => [1 - η21^2 * i for i in 0:mode1.N]))
    #     mode_op22 = DenseOperator(mode2, diagm(0 => [1 - η22^2 * i for i in 0:mode2.N]))
    #     mode_op1 = mode_op11 ⊗ mode_op12
    #     mode_op2 = mode_op21 ⊗ mode_op22
    #     Hp(t) = ion_op(t) ⊗ one(C2) ⊗ mode_op1 + one(C1) ⊗ ion_op(t) ⊗ mode_op2
    #     qoH(t) = Hp(t) + dagger(Hp(t))
    #     H = hamiltonian(T, lamb_dicke_order=0, rwa_cutoff=Inf)
    #     H1 = hamiltonian(T, lamb_dicke_order=0, rwa_cutoff=Inf, time_dependent_eta=true)
    #     H2 = hamiltonian(T, lamb_dicke_order=0, rwa_cutoff=Inf, displacement="analytic")
    #     H3 = hamiltonian(T, lamb_dicke_order=0, rwa_cutoff=Inf, displacement="analytic", time_dependent_eta=true)
    #     # only considering first order corrections to carrier (propto η^2) so this won't be perfect
    #     @test norm((qoH(tp) - H(tp, 0)).data) < 2
    #     @test norm((qoH(tp) - H1(tp, 0)).data) < 2
    #     @test norm((qoH(tp) - H2(tp, 0)).data) < 2
    #     @test norm((qoH(tp) - H3(tp, 0)).data) < 2

    #     # RWA
    #     H = hamiltonian(T, lamb_dicke_order=30, rwa_cutoff=3e5)
    #     H1 = hamiltonian(T, lamb_dicke_order=30, rwa_cutoff=3e5, time_dependent_eta=true)
    #     H2 = hamiltonian(T, lamb_dicke_order=30, rwa_cutoff=3e5, displacement="analytic")
    #     H3 = hamiltonian(T, lamb_dicke_order=30, rwa_cutoff=3e5, displacement="analytic", time_dependent_eta=true)
    #     # only considering first order corrections to carrier (propto η^2) so this won't be perfect
    #     @test norm((qoH(tp) - H(tp, 0)).data) < 2
    #     @test norm((qoH(tp) - H1(tp, 0)).data) < 2
    #     @test norm((qoH(tp) - H2(tp, 0)).data) < 2
    #     @test norm((qoH(tp) - H3(tp, 0)).data) < 2

    #     @test_throws AssertionError hamiltonian(T, lamb_dicke_order=[1, 2, 3], rwa_cutoff=3e5)
    # end
end  # end suppress
