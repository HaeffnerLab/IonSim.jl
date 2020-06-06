using QuantumOptics: displace, FockBasis, dagger, ⊗, embed, DenseOperator
using LinearAlgebra: diagm, norm
using Test, IonSim


@testset "base_functions" begin
    # get_η
    C = ca40()
    L = laser()
    chain = linearchain(
        ions=[C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[1], y=[], z=[1])
    )
    axial_mode = chain.vibrational_modes.z[1]
    radial_mode = chain.vibrational_modes.x[1]
    ηs = (2π / 7.29147e-7) * sqrt(1.0545718e-34 / (2 * 6.635943757345042e-26 * 2π))
    @test abs(get_η(axial_mode, L, C)) ≈ ηs / sqrt(2 * 1e6)
    L.k = (x=1/√2, y=0, z=1/√2)
    @test abs(get_η(axial_mode, L, C)) ≈ (ηs / sqrt(2 * 1e6)) / √2

    # _alaguerre
    L32 = (1/6) * (-2^3 + 3*(2+3)*2^2 - 3*(2+2)*(2+3)*2 + (2+1)*(2+2)*(2+3))
    @test IonSim._alaguerre(2, 3, 2) ≈ L32

    # _Dnm
    b = FockBasis(100)
    ξ = im * exp(2π*im)
    d = displace(b, ξ).data
    diff = 0.0
    for i in 1:100, j in 1:100
        diff += abs(d[i, j] - IonSim._Dnm(ξ, i, j))
    end
    @test diff < 100  # <1% difference of L1 norm
                      # Note: displace() is an approximation, whereas _Dnm should not be
    
    # _D
    d1, _ = IonSim._Ds(2, 8.9, 0.1, 9, 1, 10, 9, 1)
    @test  d1 ≈ 0.28728756457152516 - 0.49624435859832966im rtol=1e-4
    d2, d3 = IonSim._Ds(1, 0, 0.1, 1, 1, 10, 9, 1)
    @test d2 ≈ conj(d3)
    d4, d5 = IonSim._Ds(1, 0, 0.1, 1, 1, 10, 8, 1)
    @test d4 ≈ conj(d5) rtol=1e-4
end


@testset "parameter_arrays" begin
    # _Ωmatrix
    C = ca40(selected_level_structure=["S-1/2", "D-1/2"])
    C1 = ca40(selected_level_structure=["S+1/2", "D-1/2"])
    L1 = laser(); L1.pointing = [(1, 1.0), (2, 1.0)]
    L2 = laser(); L2.pointing = [(1, 1.0), (2, 1.0)]
    L1.E = 1
    L1.ϕ = 1
    L2.E = 2
    L2.ϕ = 2
    chain = linearchain(
            ions=[C, C1], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[], y=[], z=[1])
        )
    T = trap(configuration=chain, lasers=[L1, L2])
    Ωnmkj = IonSim._Ωmatrix(T, 1)
    @test 2 * real(Ωnmkj[1, 1][1](0.5)) ≈ real(Ωnmkj[1, 2][1](0.5))
    # @test 2 * real(Ωnmkj[2, 1][1](0.5)) ≈ real(Ωnmkj[2, 2][1](0.5))
    @test !( real(Ωnmkj[1, 1][1](0.5)) ≈ real(Ωnmkj[2, 1][1](0.5)) )
    
    L1.E = cos
    L1.ϕ = 0
    Ωnmkj = IonSim._Ωmatrix(T, 1)
    Ω0 = real(Ωnmkj[1, 1][1](0))
    t = 0:0.1:10
    # @test @.(real(Ωnmkj[1, 1][1](t)) / Ω0) ≈ cos.(t)

    L1.E = 1
    L1.ϕ = cos
    Ωnmkj = IonSim._Ωmatrix(T, 1)
    Ω0 = real(Ωnmkj[1, 1][1](0))
    # @test @.(Ωnmkj[1, 1][1](t)/Ω0) ≈ @.(exp(-1im * 2π * cos(t)))
    

    # _Δmatrix
    ## zero B, zero laser detuning, zero stark shift should give an array of zeros
    Δ = IonSim._Δmatrix(T, 1)
    for i in 1:2, j in 1:2
        @test Δ[i, j][1] == 0
    end
    
    ## zero B, zero stark shift, non-zero laser detuning should return columns that
    # all have the same value
    L1.Δ = 1
    L2.Δ = -1
    Δ = IonSim._Δmatrix(T, 1)
    @test Δ[1, 1][1] ≈ 2π && Δ[2, 1][1] ≈ 2π && Δ[1, 2][1] ≈ -2π && Δ[2, 2][1] ≈ -2π

    ## zero B, zero laser detuning, now just add stark shift to one of the ions
    L1.Δ = 0
    L2.Δ = 0
    C.stark_shift["S-1/2"] = 1
    Δ = IonSim._Δmatrix(T, 1)
    @test Δ[1, 1][1] ≈ -2π && Δ[1, 2][1] ≈ -2π && Δ[2, 1][1] ≈ 0 && Δ[2, 2][1] ≈ 0

    ## lastly let's test when resonant
    C.stark_shift["S-1/2"] = 0
    T.B = 1e-4
    f = transition_frequency(T, C, ("S-1/2", "D-1/2"))
    L1.Δ = f
    Δ = IonSim._Δmatrix(T, 1)
    @test Δ[1, 1][1] ≈ 0


    # _ηmatrix
    chain = linearchain(
                ions=[C, C1], com_frequencies=(x=2e6,y=2e6,z=1e6), selected_modes=(x=[1], y=[2], z=[1])
            )
    L1.k = (x=1/√2, y=0, z=1/√2)
    L2.k = (x=0, y=1/√2, z=1/√2)
    L3 = laser(); L3.pointing = [(1, 1.0), (2, 1.0)]
    L3.k = (x=0, y=0, z=1)
    T = trap(configuration=chain, lasers=[L1, L2, L3])
    η = IonSim._ηmatrix(T, 1e-6)
    @test abs(η[1, 1, 1](1.0)) ≈ abs(η[1, 1, 3](1.0) / √2)
    @test η[1, 1, 1](1.0) ≈ η[2, 1, 1](1.0)
    @test η[1, 2, 2](1.0) ≈ -η[2, 2, 2](1.0)
    @test η[1, 3, 1](1.0) ≈ 0 && η[1, 3, 2](1.0) ≈ 0 && η[1, 3, 3](1.0) !== 0
    chain.vibrational_modes.z[1].δν = t -> 1e6t
    η = IonSim._ηmatrix(T, 1e-6)
    @test η[1, 1, 3](0.0) ≈ 2 * η[1, 1, 3](3) 
end

@testset "setup_functions" begin
    # _setup_global_B_hamiltonian
    C = ca40()
    L = laser()
    chain = linearchain(
        ions=[C], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[], y=[], z=[1])
    )
    T = trap(configuration=chain, lasers=[L], δB=0)
    global_B_indices, global_B_scales, bfunc = IonSim._setup_global_B_hamiltonian(T, 1)

    ## T.δB = 0 -> global_B_indices and global_B_scales should be empty arrays 
    @test length(global_B_indices) == 0 && length(global_B_scales) == 0

    T.δB = sin 
    t = 0:0.1:10
    T.configuration.vibrational_modes.z[1].N = 3
    global_B_indices, global_B_scales, bfunc = IonSim._setup_global_B_hamiltonian(T, 1)
    indxs = [[i + 8j for j in 0:2] for i in 1:8]
    matched = unique([indxs; global_B_indices])
    @test length(matched) == 8
    @test bfunc.(t) == 2π .* sin.(t)
    @test global_B_scales[end] == zeeman_shift(1, C.selected_level_structure["D+5/2"]) 
    

    # _setup_δν_hamiltonian
    ## should return empty arrays if δν=0
    δν_indices, δν_functions = IonSim._setup_δν_hamiltonian(T, 1)
    @test length(δν_indices) == 0 && length(δν_functions) == 0
    
    T.configuration.vibrational_modes.z[1].N = 10
    T.configuration.vibrational_modes.z[1].δν = 1
    indxs = [[8j + i for i in 1:8] for j in 1:9]
    δν_indices, δν_functions = IonSim._setup_δν_hamiltonian(T, 1)
    matched = unique([indxs; δν_indices])
    @test length(matched) == 10
    @test length(δν_functions) == 1

    chain = linearchain(
        ions=[C], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[1], y=[], z=[1])
    )
    T = trap(configuration=chain, lasers=[L], δB=0)
    T.configuration.vibrational_modes.x[1].δν = cos
    T.configuration.vibrational_modes.z[1].δν = sin
    δν_indices, δν_functions = IonSim._setup_δν_hamiltonian(T, 1)
    t = 0:0.1:10
    @test δν_functions[1].(t) == 2π .* cos.(t)
    @test δν_functions[2].(t) == 2π .* sin.(t)


    # _setup_fluctuation_hamiltonian
    T = trap(configuration=chain, lasers=[L], δB=1)
    T.configuration.vibrational_modes.x[1].δν = cos
    T.configuration.vibrational_modes.z[1].δν = sin
    δν_indices, δν_functions = IonSim._setup_δν_hamiltonian(T, 1)
    global_B_indices, global_B_scales, bfunc = IonSim._setup_global_B_hamiltonian(T, 1)
    all_unique_indices, gbi, gbs, bfunc1, δνi, δνfuncs = IonSim._setup_fluctuation_hamiltonian(T, 1)
    @test any(gbi .== global_B_indices)
    @test any(gbs .== global_B_scales)
    @test bfunc(17.0) == bfunc1(17.0)
    @test any(δνi .== δν_indices)
    @test δνfuncs[end](92.0) == δν_functions[end](92.0)


    # _setup_base_hamiltonian
    ## we can check that the functional description is correct when test hamiltonian and/or
    ## the time-evolution. For now let' just make sure elements are getting assinged to the 
    ## correct places. We'll look at two 2-level ions coupled to a single mode with dim=2.
    ## Note: qo embed performs tensor product in backwards order compared to what we'd usually
    ## expect
    C = ca40(selected_level_structure=["S-1/2", "D-1/2"])
    L = laser()
    chain = linearchain(ions=[C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[], y=[], z=[1]))
    T = trap(configuration=chain, B=4e-4, Bhat=(x=0,y=0,z=1), lasers=[L])
    mode = T.configuration.vibrational_modes.z[1]
    mode.N = 2

    ## first just shine light on 1st ion 
    L.pointing = [(1, 1.0)]
    _, repeated_indices, conj_repeated_indices = IonSim._setup_base_hamiltonian(T, 1, 100, 1e10)
    ion1_indxs = [(2,1), (4,3), (2,5), (4,7), (6,1), (8,3), (8,7), (6,5)]
    unique_els = unique([IonSim._flattenall(repeated_indices); ion1_indxs])
    # @test length(repeated_indices) == 4
    # @test length(unique_els) == 8

    ## now just shine light on second ion
    L.pointing = [(2, 1.0)]
    _, repeated_indices, conj_repeated_indices = IonSim._setup_base_hamiltonian(T, 1, 100, 1e10)
    ion2_indxs = [(3,1), (4,2), (4,6), (3,5), (7,1), (8,2), (7,5), (8,6)]
    unique_els = unique([IonSim._flattenall(repeated_indices); ion2_indxs])
    # @test length(repeated_indices) == 4
    # @test length(unique_els) == 8

    ## now shine light on both ions
    L.pointing = [(1, 1.0), (2, 1.0)]
    _, repeated_indices, conj_repeated_indices = IonSim._setup_base_hamiltonian(T, 1, 100, 1e10)
    unique_els = unique([IonSim._flattenall(repeated_indices); ion1_indxs; ion2_indxs])
    # @test length(repeated_indices) == 8
    # @test length(unique_els) == 16

    ## now lets set lamb_dicke_order = 0
    _, repeated_indices, conj_repeated_indices = IonSim._setup_base_hamiltonian(T, 1, 0, 1e10)
    ld_indxs = [(2,1), (4,3), (3,1), (4,2), (8,7), (6,5), (7,5), (8,6)]
    unique_els = unique([IonSim._flattenall(repeated_indices); ld_indxs])
    # @test length(repeated_indices) == 4
    # @test length(unique_els) == 8

    ## when laser tuned to carrier, setting an rwa_cutoff below the vibrationl_mode frequency should hvae the same affect
    L.Δ = transition_frequency(T, C, ("S-1/2", "D-1/2"))
    _, repeated_indices, conj_repeated_indices = IonSim._setup_base_hamiltonian(T, 1, 100, 1e5)
    unique_els = unique([IonSim._flattenall(repeated_indices); ld_indxs])
    # @test length(repeated_indices) == 4
    # @test length(unique_els) == 8

    ## finally lets make sure the conj_repeated_indices are going to the right places
    _, repeated_indices, conj_repeated_indices = IonSim._setup_base_hamiltonian(T, 1, 100, Inf)
    expected_conj_repeated_indices = [(2,5), (4,7), (4,6), (3,5)]
    unique_els = unique([IonSim._flattenall(conj_repeated_indices); expected_conj_repeated_indices])
    # @test length(repeated_indices) == 6
    # @test length(unique_els) == 5  # 4 actual indices plus (-1,0) flag
end


@testset "hamiltonian" begin
    C = ca40(selected_level_structure=["S-1/2", "D-1/2"])
    L = laser()
    L.pointing = [(1, 1.0), (2, 1.0)]
    chain = linearchain(ions=[C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), selected_modes=(x=[], y=[], z=[1]))
    T = trap(configuration=chain, B=4e-4, Bhat=(x=0,y=0,z=1), lasers=[L])
    mode = T.configuration.vibrational_modes.z[1]
    L.Δ = transition_frequency(T, C, ("S-1/2", "D-1/2")) + 2e5
    L.ϕ = 1/2
    mode.N = 30
    Efield_from_rabi_frequency!(0.5e6, T, 1, 1, ("S-1/2", "D-1/2"))

    
    # full hamiltonian (w/o conj_repeated_indices)
    # full hamiltonian (w/o conj_repeated_indices)
    ion1_op(t) = π/2 * exp(-im * 2π *(0.2t + 1/2)) * C["D-1/2"] ⊗ C["S-1/2"]'
    ion2_op(t) = π/2 * exp(-im * 2π * (0.2t + 1/2)) * C["D-1/2"] ⊗ C["S-1/2"]'
    η = IonSim.get_η(mode, L, C)
    mode_op(t) = displace(mode.basis, im * η * exp(im * 2π * t))
    Hp(t) = embed(get_basis(T), [1, 3], [ion1_op(t), mode_op(t)]) + embed(get_basis(T), [2, 3], [ion2_op(t), mode_op(t)])
    qoH(t) = Hp(t) + dagger(Hp(t))
    H = hamiltonian(T, lamb_dicke_order=30, rwa_cutoff=1e10)
    # @test norm(real.((qoH(1027.32) - H(1027.32, 0)).data[1:10, 1:10])) < 1e-6
    # @test norm(imag.((qoH(1027.32) - H(1027.32, 0)).data[1:10, 1:10])) < 1e-6


    # full hamiltonian (w/ conj_repeated_indices)
    H = hamiltonian(T, lamb_dicke_order=30, rwa_cutoff=Inf)
    # @test norm(real.((qoH(1027.32) - H(1027.32, 0)).data[1:10, 1:10])) < 1e-6
    # @test norm(imag.((qoH(1027.32) - H(1027.32, 0)).data[1:10, 1:10])) < 1e-6
    
    # Lamb-Dicke
    mode_op(t) = DenseOperator(mode.basis, diagm(0 => [1 - η^2 * i for i in 0:29]))
    Hp(t) = embed(get_basis(T), [1, 3], [ion1_op(t), mode_op(t)]) + embed(get_basis(T), [2, 3], [ion2_op(t), mode_op(t)])
    qoH(t) = Hp(t) + dagger(Hp(t))
    H = hamiltonian(T, lamb_dicke_order=0, rwa_cutoff=Inf)
    # only considering first order corrections to carrier (propto η^2) so this won't be perfect
    # @test norm(real.((qoH(1027.32) - H(1027.32, 0)).data[1:10, 1:10])) < 0.05  
    # @test norm(imag.((qoH(1027.32) - H(1027.32, 0)).data[1:10, 1:10])) < 0.05

    # RWA
    H = hamiltonian(T, lamb_dicke_order=30, rwa_cutoff=3e5)
    # only considering first order corrections to carrier (propto η^2) so this won't be perfect
    # @test norm(real.((qoH(1027.32) - H(1027.32, 0)).data[1:10, 1:10])) < 0.05  
    # @test norm(imag.((qoH(1027.32) - H(1027.32, 0)).data[1:10, 1:10])) < 0.05
end