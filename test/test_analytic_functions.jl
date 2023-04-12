using Test, IonSim, IonSim.analytical

@testset "analytic_functions -- MolmerSorensen" begin
    t = 0:100:10000
    ν = 1
    Ω = 0.025
    δ = 0.95
    η = 0.1
    n̄ = 0

    # test off-resonant limit: ηΩ << (ν - δ)
    # should get rabi oscillations between |ee> and |gg> with Ω̃ = (ηΩ)^2 / (ν - δ)
    ρgg, ρee = molmersorensen2ion(t, Ω, ν, δ, η, n̄)
    @test ρee ≈ @.(sin(2π * (η * Ω)^2 / 2(ν - δ) * t)^2)

    # test near resonant limit
    # if we set Ω = (ν-δ) / (2η√K); K<:Int, we should get ρee=ρgg=0.5 at time τ = √K / (2ηΩ)
    Ω = (ν - δ) / (2η)
    τ = 1 / (2η * Ω)
    ρgg, ρee = molmersorensen2ion([τ], Ω, ν, δ, η, n̄)
    @test ρgg[1] + ρee[1] ≈ 1.0

    Ω = (ν - δ) / (2η * √2)
    τ = √2 / (2η * Ω)
    ρgg, ρee = molmersorensen2ion([τ], Ω, ν, δ, η, n̄)
    @test ρgg[1] + ρee[1] ≈ 1.0

    # test hot gate
    Ω = (ν - δ) / (2η)
    τ = 1 / (2η * Ω)
    n̄ = 50
    ρgg, ρee = molmersorensen2ion([τ / 2, τ], Ω, ν, δ, η, n̄)
    @test ρgg[2] + ρee[2] ≈ 1.0 rtol = 1e-2
    @test ρgg[1] + ρee[1] ≈ 0.75 rtol = 1e-2
end

@testset "analyticfunctions -- RabiFlop" begin
    # # test carrier with motion in ground state
    t = 0:0.01:1
    Ω = 1
    n̄ = 0
    η = 0
    @test rabiflop(t, Ω, η, n̄) ≈ @.(sin(2π * t / 2)^2) rtol = 1e-5

    # test bluesideband, Ω̃ ≈ ηΩ
    η = 0.1
    @test rabiflop(t, Ω, η, n̄, s=1) ≈ @.(sin(2π * η * t / 2)^2) rtol = 1e-2
end
