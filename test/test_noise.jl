using Test, IonSim, IonSim.noise
using FFTW
using StochasticDiffEq

@testset "noise -- GaussianProcess" begin

    # == Exact Power Spectrums for ARMA process == #
    psd_ar1 = (f, ϕ, σ, dt) -> σ^2 / (1 + ϕ^2 - 2 * ϕ * cos(2pi * f * dt))
    # ar_psd(sigma, phi, freq, T) = sigma^2 ./ map(f -> abs.(1.0 - dot(phi, exp.(-1im * (2pi * f * T) .* (1:length(phi))))).^2, freq)

    dt = 1
    tscale = 1e-6
    ϕ = 0.5
    OU = GaussianProcess(f -> psd_ar1(f, ϕ, 1, dt * tscale), dt, tscale; coeff_length = 1)
    # OU = GaussianProcess(f -> ar_psd(1.0, ϕ, f, dt*tscale), dt, tscale; coeff_length=2)
    arma = OU.ARMAProcess.arma

    # == Power Spectrum == #
    OU_power = map(f -> psd_ar1(f, ϕ, arma.sigma, dt * tscale), arma.freq) #ar_psd(arma.sigma, arma.phi, arma.freq, arma.dt*arma.timescale)
    arma_power = arma.psd

    # plot(arma.freq, OU_power, label="OU Model")
    # plot!(arma.freq, arma_power, label="GaussianProcess Fit")

    # OU Parameters
    θ = (-log(ϕ) / dt)
    μ = 0
    σ = (arma.sigma / √(dt)) * (2 * log(ϕ) / (ϕ^2 - 1))^(1 / 2)

    dt = 1
    tspan = 0:dt:(dt * 1000)
    est_avg = zeros(length(tspan))
    exact_avg = zero(est_avg)
    Ntrials = 1024

    for i in 1:Ntrials
        OU_est = OU(0, 0, continuous = true) # returns W(t)
        _, w = calculate_noise!(OU_est, tspan) # return dW/dt = x(t)
        est_avg += (dt / (tspan[end] / dt)) .* abs.(fft(w)) .^ 2 ./ Ntrials

        OU_exact = StochasticDiffEq.OrnsteinUhlenbeckProcess(θ, μ, σ, 0.0, 0.0)
        OU_exact.dt = dt
        map(t -> OU_exact(t), tspan) # returns x(t)
        exact_avg += (dt / (tspan[end] / dt)) .* abs.(fft(OU_exact.u)) .^ 2 ./ Ntrials
        # plot!(t,w)
    end

    N = length(est_avg) ÷ 2
    freq = fftfreq(2N, 1 / (dt * tscale))
    # plot!(freq[1:N], est_avg[1:N], label="GaussianProcess (CARMA) Ensemble Avg")
    # plot!(freq[1:N], exact_avg[1:N], label="OrnsteinUhlenbeckProcess Ensemble Avg")
    # plot!(xlimit=(0, min(maximum(arma.freq), maximum(freq)))) |> display

    # # == ACF == #
    time = 0:20
    # plot(time, real.(ifft(OU_power))[1:length(time)],linewidth=2)
    # plot!(time, real.(ifft(arma_power))[1:length(time)]) |> display

    @test isapprox(OU_power, arma_power, rtol = 1e-10)
    @test isapprox(
        real.(ifft(OU_power))[1:length(time)],
        real.(ifft(arma_power))[1:length(time)],
        rtol = 1e-10
    )
    @test isapprox(est_avg, exact_avg, rtol = 0.1)
end
