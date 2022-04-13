using FFTW
using LinearAlgebra
using IterativeSolvers
using PyCall
np = pyimport("numpy")
so = pyimport("scipy.optimize")

# == AR Power Spectrum Estimate == #
ar_psd(sigma, phi, freq, T) = sigma^2 ./ map(f -> abs.(1.0 - dot(phi, exp.(-1im * (2pi * f * T) .* (1:length(phi))))).^2, freq)

mutable struct ARMA
    # Paramters
    psd::Vector      # (Discrete) Power Spectral Density
    ar_psd::Vector   # AR Power Spectral Denisty Estimate
    freq::Union{Vector,Frequencies} # Frequency values
    acf::Vector      # Ensemble Autocorrelation Function
    ar_acf::Vector   # AR ACF Estimate
    ma_acf::Vector   # MA ACF Estimate
    phi::Vector      # AR parameters phi_1, ..., phi_p
    psi::Vector      # MA parameters theta_1, ..., theta_q
    ar_roots::Vector # AR roots for each phi_i

    # Extra
    mean::Real       # Mean of noise process
    sigma::Real      # STD estimate of white noise
    dt::Real         # Correlation timestep
    normalize_sigma::Bool # process variance normalized
    coeff_length::Integer # AR coefficient length
    impulse_length::Integer # MA Coefficient length
    freq_bin::Real   # Frequency bin size (cycle / (window_length * dt))
end

function ARMA(psd::AbstractVector, freq::Union{AbstractVector,Frequencies}; mean::Real=0.0, normalize_sigma::Bool=true,  
    coeff_length::Union{Integer,Nothing}=nothing, impulse_length::Union{Integer,Nothing}=nothing, dt::Real=1.0, timescale::Real=1.0)
       
    # == Autocorrelation Function == #
    w_triangle = 1.0 .- (0:length(psd)-1)/length(psd)
    acf = (1 / dt) .* real.(ifft(psd)) .* w_triangle # Normalize to bin size = N * dt

    # == Find ARMA Coefficients == #
    if coeff_length === nothing
      coeff_length = length(acf)
    end
    if impulse_length === nothing
      impulse_length = length(acf)
    end
    phi, psi, ar_roots = generate_ARMA_coefficients(acf, coeff_length, impulse_length)
    
    # == Noise amplitude estimate == #
    sigma = sqrt(acf[1] - phi' * acf[2:length(phi)+1])
    
    # == AR Fitted PSD & ACF == #
    # sigma, _ = so.curve_fit((f,σ...)->ar_psd(σ[1], phi,f,dt), abs.(freq), psd, sigma)
    # sigma = sigma[1]

    ar_ps = ar_psd(sigma, phi, abs.(freq), dt * timescale)
    ar_acf = real.(ifft(ar_ps)) .* w_triangle
    ma_acf = sigma^2 .* [dot(psi[1:end-i], psi[1+i:end]) for i in 0:(length(psi)-1)]

    # == Update noise amplitude est. & normalize == #
    if normalize_sigma # process with unit variance
      sigma /= sqrt(ar_acf[1])
      ar_ps ./= ar_acf[1]
      ma_acf ./= ar_acf[1]
      ar_acf ./= ar_acf[1]
      psd ./= acf[1]
      acf ./= acf[1]
    end


    return ARMA(psd, ar_ps, freq, acf, ar_acf, ma_acf, phi, psi, ar_roots, mean, sigma, dt, normalize_sigma, coeff_length, impulse_length, abs(freq[2] - freq[1]))
end

# Compute the impulse response function associated with ARMA process arma
function generate_IIR_coefficients(phi::Vector, impulse_length)
  p = length(phi)

  if impulse_length === nothing
      impulse_length = p
  else
      err_msg = "Impulse length must be greater than number of AR coefficients"
      @assert impulse_length >= p err_msg
  end

  # == Pad theta with zeros at the end == #
  psi_zero = 1.0
  psi = zeros(impulse_length)
  for j = 1:(impulse_length)
      for i = 1:min(j, p)
          psi[j] += phi[i] * (j-i > 0 ? psi[j-i] : psi_zero)
      end
  end

  return [psi_zero; psi[1:end-1]]
end

function generate_AR_coefficients(acf::Vector, coeff_length)
  
  if !(coeff_length <= length(acf) - 1)
    @assert false "coeff_length must be <= (length(acf) - 1) in Yule-Walker equation"
  end

  # == AR Coefficients Solver == # 
  ρ = acf[1:coeff_length+1] / acf[1]
  n = length(ρ) - 1
  Rmat = zeros(n,coeff_length)
  for j in 1:n
    Rmat[j,:] .= ρ[abs.(j .- (1:coeff_length)) .+ 1]
  end
  phi = IterativeSolvers.lsqr(Rmat,ρ[2:end])

  charpoly = [-reverse(phi);1.0]
  roots = np.roots(charpoly) # roots of |z|
  if any(abs.(roots) .<= 1.0) # |z| > 1 s.t. |w| < 1 to be stationary
    @warn "AR(p) process is not stationary, must satisfy characteristic equation with |z| > 1. Consider trying different coefficient length."
  end

  return phi, roots
end

  # for j in 1:coeff_length
  #     Rmat[j,j] = 1.0
  #     if j < coeff_length
  #       Rmat[j,j+1:end] = ρ[1:end-j]
  #     end
  # end

function generate_ARMA_coefficients(acf::Vector, coeff_length, impulse_length)
  # == AR Coefficients (FIR Filter) == #
  phi, ar_roots = generate_AR_coefficients(acf, coeff_length)

  # == MA Coefficients (IIR Filter) == #
  psi = generate_IIR_coefficients(phi, impulse_length)

  return (phi, psi, ar_roots)
end
