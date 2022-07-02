using FFTW
using PyCall
using Plots
using QuantumOptics
using StochasticDiffEq
np = pyimport("numpy")
include("GaussianNoiseProcess.jl")

# == Desired PSD Function == #
function color_spectrum(f, f0, α, dt)
    # fnyquist = 1/(2dt)
    return 1 / (2pi * (abs(f) < f0 ? f0 : abs(f)) * dt)^α
end

# == Generate Noise Process == #
tscale = 1e-3 #ms
dt = 2 * tscale #0.002ms steps (2μs)
f0 = 2.0e3
α = 1
PinkNoise = GaussianNoiseProcess.GaussianProcess(
    f -> color_spectrum(f, f0, α, dt * tscale),
    dt,
    tscale;
    window_length = 2^10,
    coeff_length = 2^6
)
arma = PinkNoise.ARMAProcess.arma; # for simplicity, contains all ARMA info

# == Test noise == #
Δt = dt # if you try to sample at Δt != dt, then switch to using PinkNoise(0,0,method="CARMA")
tspan = 0:Δt:2.0
N = length(tspan)

Ntrials = 512
psd_arma = map(
    i -> begin
        W = PinkNoise(0, 0) # Creates NoiseProcess that intercaces with solver
        x = GaussianNoiseProcess.calculate_noise!(W, tspan)[2] # interface that returns noise
        (W.dt / dt) * (1 / (tspan[end] / W.dt)) .* abs.(np.fft.fft(x)) .^ 2 ./ Ntrials # Power Spectrum Est.
    end,
    1:Ntrials
)

psd_arma = sum(hcat(psd_arma...), dims = 2)[:] # averaging over power spectrum
freq = fftfreq(N, 1 / ((tspan[2] - tspan[1]) * tscale))

# == Plotting PSD == #
plot(freq[2:end], psd_arma[2:end], label = "Ensemble Avg")
plot!(arma.freq[2:end], arma.psd[2:end], label = "Model")
plot!(arma.freq[2:end], arma.arma_psd[2:end], label = "ARMA Est.")
plot!(xscale = :log, yscale = :log)
plot!(xlimit = (1e2, 1e6), ylimit = (1e-3, Inf)) |> display

# == Plotting ACF == #
T = 0:1:(length(arma.acf) - 1)
T *= dt
plot(T, arma.acf, label = "Model", linewidth = 2.0)
plot!(T, arma.arma_acf, label = "ARMA Est.")

ens_acf = (dt / Δt) * np.fft.irfft(psd_arma, length(psd_arma)) # Julia FFTW doesn't support this function
L = length(ens_acf) ÷ 2 + 1
T = 0:1:(L - 1)
T *= Δt
plot!(T, ens_acf[1:L], linewidth = 2, label = "Ensemble Avg")
# plot!(xlimit=(0,1))
plot!(xlabel = "time (ms)") |> display

# == Stochastic TISE == #

# construct bases
sb = SpinBasis(1 // 2)
vb = FockBasis(10)
cb = sb ⊗ vb;

# Electronic ground state in Z-basis
ψg = spindown(sb)
# Electronic excited state in Z-basis
ψe = spinup(sb)
# Inital state
ψi = ψg ⊗ fockstate(vb, 0);

# Bare rotated-basis, single-ion, VAET Hamiltonian
H(J, Δ, κ, ν) =
    2π * (J / 2) * sigmax(sb) ⊗ one(vb) +
    2π * (Δ / 2) * sigmaz(sb) ⊗ one(vb) +
    2π * (κ / 2) * sigmaz(sb) ⊗ (create(vb) + destroy(vb)) +
    2π * ν * one(sb) ⊗ number(vb)

# Measurement function
fout(t, ψ) = real.(expect((ψe ⊗ ψe') ⊗ one(vb), ψ))

# σz noise operator for use with Lindblad master equation
J_Δ(Γ) = [√(Γ) * sigmaz(sb) ⊗ one(vb)]
J_ν(Γ) = [√(Γ) * one(sb) ⊗ number(vb)]

# setup
tspan = 0:0.01:2.0
Δt = tspan[2] - tspan[1]

# notice that I chose a timescale for the noise and the frequencies
params = (J = 1.3, Δ = -1.2, κ = 0.229, ν = 1.8) #kHz
Γ = 1.0 # arb.

Ntrials = 256 # minimum number to get to ~10-20% error
alg = StochasticDiffEq.RKMil() # found to be faster and more accurate
adaptive = false # should be the default
dt = 1e-4 # tested to be the minimum step needed to get accurate results out to 2ms

# Quantum trajectories average
sol =
    map(
        trial -> stochastic.schroedinger(
            tspan,
            ψi,
            H(params.J, params.Δ, params.κ, ν),
            J_Δ(Γ);
            fout = fout,
            noise = PinkNoise(0, 0),
            dt = dt,
            adaptive = adaptive,
            alg = alg,
            reltol = 1e-8,
            abstol = 1e-8,
            normalize_state = true
        )[2],
        1:Ntrials
    ) |> vec -> hcat(vec...) |> mat -> sum(mat, dims = 2) ./ Ntrials

plot(tspan, sol)
