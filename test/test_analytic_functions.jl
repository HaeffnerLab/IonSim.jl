using Test, IonSim, IonSim.analytical


@testset "MolmerSorensen" begin
    t = 0:100:10000
    ν = 1
    Ω = 0.025
    δ = 0.95
    η = 0.1
    n̄ = 0

    # test off-resonant limit: ηΩ << (ν - δ)
    # should get rabi oscillations between |ee> and |gg> with Ω̃ = (ηΩ)^2 / (ν - δ)
    ρgg, ρee = two_ion_ms(t, Ω, ν, δ, η, n̄)
    @test ρee ≈ @.(sin(2π * (η * Ω)^2 / 2(ν - δ) * t)^2)
    
    # test near resonant limit
    # if we set Ω = (ν-δ) / (2η√K); K<:Int, we should get ρee=ρgg=0.5 at time τ = √K / (2ηΩ)
    Ω = (ν - δ) / (2η)
    τ = 1 / (2η * Ω)
    ρgg, ρee = two_ion_ms([τ], Ω, ν, δ, η, n̄)
    @test ρgg[1] + ρee[1] ≈ 1.0

    Ω = (ν - δ) / (2η * √2)
    τ = √2 / (2η * Ω)
    ρgg, ρee = two_ion_ms([τ], Ω, ν, δ, η, n̄)
    @test ρgg[1] + ρee[1] ≈ 1.0

    # test hot gate
    Ω = (ν - δ) / (2η)
    τ = 1 / (2η * Ω)
    n̄ = 50
    ρgg, ρee = two_ion_ms([τ/2, τ], Ω, ν, δ, η, n̄)
    @test ρgg[2] + ρee[2] ≈ 1.0 rtol=1e-2
    @test ρgg[1] + ρee[1] ≈ 0.75 rtol=1e-2
end


@testset "RabiFlop" begin
    # TODO: uptdate tests for new rabi_flop function 

    # # test carrier with motion in ground state
    # t = 0:0.01:1
    # Ω = 1
    # n̄ = 0
    # η = 0.1
    # @test rabi_flop(t, Ω , η, n̄) ≈ @.(sin(2π * t)^2)

    # # test bluesideband, Ω̃ ≈ ηΩ
    # @test rabi_flop(t, Ω , η, n̄, blue_sideband=true) ≈ @.(sin(2π * η * t)^2) rtol=5e-2

    # # test hot carrier, η²n̄ = 0.2 should reduce maximum excitation by ≈ 0.1
    # n̄ = 20
    # @test maximum(rabi_flop(t, Ω , η, n̄)) ≈ 0.9 rtol=1e-3

    # # test multiple vibrational modes 
    # @test rabi_flop(t, 1, η, 2n̄) ≈ rabi_flop(t, 1, [η, η], [n̄, n̄])

end
