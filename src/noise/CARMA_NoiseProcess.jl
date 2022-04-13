using FFTW
using DiffEqNoiseProcess
using LinearAlgebra: dot
using Random
using PyCall
so = pyimport("scipy.optimize")

@inline wiener_randn(rng::AbstractRNG,::Type{T}) where T = randn(rng,T)

# == ARMAProcess: Generates ARMA Coefficients and Sinc-interpolated NoiseProcess == #

arma_psd(arma, f, T) = arma.sigma^2 / abs.( 1 - dot(arma.phi, exp.(-1im * (2pi * f * T) .* (1:1:length(arma.phi))))).^2

@inline function cont_autocorr(τ,α,κ)
  α = α[1:length(α)÷2] .+ 1im.*α[length(α)÷2+1:end]
  res = (α * transpose(α)) # d[i] * d[j] (row = i, column = j)
  res ./= (κ .+  transpose(κ)) # κ[i] + κ[j]
  res = -exp.(τ * transpose(κ)) * res
  return sum(res, dims=2)[:]
end

@inline function jac_autocorr(τ,α,κ,σ)
  α = α[1:length(α)÷2] .+ 1im.*α[length(α)÷2+1:end]

  Γik = α ./ (transpose(κ) .+ κ)
  exp_λτ = exp.(transpose(κ) .* τ)

  dγdc1 = exp_λτ .* sum(Γik,dims=1)
  dγdc1 .+= exp_λτ * Γik

  dγdc2 = 1im .* dγdc1

  dγdc1 = [real.(dγdc1); imag.(dγdc1)]
  dγdc2 = [real.(dγdc2); imag.(dγdc2)]
  
  return -σ^2 .* hcat(dγdc1,dγdc2)
end

# @inline function cont_autocorr(τ,d,κ)
#   d = d[1:length(d)÷2] .+ 1im.*d[length(d)÷2+1:end]
#   res = (d' .* d) # d[i] * d[j] (row = i, column = j)
#   res ./= (κ' .+  κ) # κ[i] + κ[j]
#   res = -exp.(τ * transpose(κ)) * res
#   return sum(res, dims=2)[:]
# end

function ReImCoeffSolver(x, α, κ, σ)
  t = x[1:length(x)÷2]

  # == Re and Im solution == #
  res = @fastmath(cont_autocorr(t, α, κ))
  
  return σ^2 .* [real.(res); imag.(res)]
end

struct CARMA
  κ::Vector
  α::Vector
  acf::Vector
end

struct CARMAProcess
  psd_model::Function
  dt::Real
  Hfreq::Function
  Htime::Function
  normalize_sigma::Bool
  carma::CARMA
  arma::ARMA
  ARMAProcess::ARMAProcess
  window_length::Int
  timescale::Real

  @inline function CARMAProcess(psd::Function, dt::Real; mean::Real=0.0, normalize_sigma::Bool=true, Hfreq::Union{Function,Nothing}=nothing, Htime::Union{Function,Nothing}=nothing, 
    window_length::Integer=2^10, coeff_length::Union{Integer,Nothing}=2^6, impulse_length::Union{Integer,Nothing}=nothing, timescale::Real=1.0,
    kwargs...)
    
    discrete = ARMAProcess(psd, dt; mean, normalize_sigma, Hfreq, Htime, 
      window_length, coeff_length, impulse_length, timescale,
      kwargs...)

    if all(isapprox.(1.0 .+ discrete.arma.phi, 1.0))
      carma_κ = carma_α = carma_acf = []
    else
      # == continous model poles == #
      arma_ρ = abs.(discrete.arma.ar_roots)
      arma_θ = atan.(imag.(discrete.arma.ar_roots) ./ real.(discrete.arma.ar_roots))
      carma_κ = -(log.(arma_ρ) .+ 1im .* arma_θ)
        
      # == Fit Model == #
      ts = 0:1:(2*length(discrete.arma.phi)-1)
      ys = discrete.arma.acf[ts .+ 1]

      # acf_iir = discrete.arma.sigma^2 .* [dot(discrete.arma.psi[1:end-i], discrete.arma.psi[1+i:end]) for i in 0:(length(discrete.arma.psi)-1)]
      # ys = discrete.arma.ma_acf[ts .+ 1]
      carma_α, _ = so.curve_fit((x,d...)->ReImCoeffSolver(x,collect(d),carma_κ, discrete.arma.sigma), [ts;ts], [ys;zeros(length(ys))], 
                    ones(Float64, 2*length(discrete.arma.phi)), jac=(x,d...)->jac_autocorr(x[1:length(x)÷2],collect(d),carma_κ,discrete.arma.sigma));
      #! add error statement if doesn't fit good
      ts = 0:(length(discrete.arma.acf)-1)
      carma_acf = ReImCoeffSolver(ts, carma_α, carma_κ, discrete.arma.sigma)
      carma_α = carma_α[1:length(discrete.arma.phi)] .+ 1im .* carma_α[length(discrete.arma.phi)+1:end]
    end

    # == Return ARMAProcess structure == #
    new(discrete.psd, discrete.dt, discrete.Hfreq, discrete.Htime, discrete.normalize_sigma, CARMA(carma_κ, carma_α, carma_acf), discrete.arma, discrete, discrete.window_length, discrete.timescale)

  end
end

# == Generators of Discrete AR/MA time series == #
struct CARMA_Process{T1,T2,T3,T4,T5,T6,T7,T8}
  κ::T1
  α::T2
  sigma::T3
  mean::T4
  dt::T5
  Xpast::T6
  p::T7
  filter::T8
end

# == Helper functions == #
# et(rng,T) = wiener_randn(rng, T) + 1im * wiener_randn(rng, T)
@inline step(mean, std, rng, T) = mean .+ std .* (wiener_randn(rng, T) + 1im * wiener_randn(rng, T))

@inline function (X::CARMA_Process)(dW,W,dt,u,p,t,rng)
  
  # == calculate step == #
  T = (typeof(dW) <: AbstractArray) ? dW : typeof(dW) #? Still necessary?

  mean = @fastmath X.Xpast .* exp.(X.κ * dt)
  std = @fastmath X.sigma .* sqrt.(1.0 .- exp.(2.0 .* X.κ .* dt))
  # std = X.sigma .* X.α .* sqrt.((1.0 .- exp.(2.0 .* X.κ .* dt)) ./ (-2.0 .* X.κ))
  # std = X.sigma .* sqrt(dt / X.dt) .* sqrt.(1.0 .- exp.(2.0 .* X.κ .* X.dt))

  Wr = step(mean, std, rng, T)
  Wf = sum(real.(Wr))
  while X.filter(Wf)
    Wr = step(mean, std, rng, T)
    Wf = sum(real.(Wr))
  end
  # setproperty!(X,:Xpast, Wr)
  # Wf = sum(real.(Wr .- X.Xpast))
  X.Xpast .= Wr

  return Wf * dt
  # return Wf * dt - W[end]
  # return Wf * dt
  # return sum(real.(X.κ .* Wr)) * dt
end

# @inline function (X::CARMA_Process)(dW,W,dt,u,p,t,rng)
  
#   # == calculate step == #
#   T = (typeof(dW) <: AbstractArray) ? dW : typeof(dW) #? Still necessary?

#   mean = @fastmath X.Xpast .* exp.(X.κ * dt)
#   std = @fastmath X.sigma .* sqrt.(1.0 .- exp.(2.0 .* X.κ .* dt))
#   # std = X.sigma .* X.α .* sqrt.((1.0 .- exp.(2.0 .* X.κ .* dt)) ./ (-2.0 .* X.κ))
#   # std = X.sigma .* sqrt(dt / X.dt) .* sqrt.(1.0 .- exp.(2.0 .* X.κ .* X.dt))

#   Wr = step(mean, std, rng, T)
#   Wf = sum(real.(Wr))
#   while X.filter(Wf)
#     Wr = step(mean, std, rng, T)
#     Wf = sum(real.(Wr))
#   end
#   # setproperty!(X,:Xpast, Wr)
#   X.Xpast .= Wr

#   # return Wf - W[end]
#   return Wf * dt - W[end]
#   # return Wf * dt
# end

# @inline function (X::CARMA_Process)(dW,W,dt,u,p,t,rng)
  
#   # == calculate step == #
#   T = (typeof(dW) <: AbstractArray) ? dW : typeof(dW) #? Still necessary?

#   mean = @fastmath X.Xpast .* (exp.(X.κ * dt) .- 1.0)
#   std = @fastmath X.sigma .* sqrt.(1.0 .- exp.(2.0 .* X.κ .* dt))

#   dWr = step(mean, std, rng, T)
#   dWf = sum(real.(dWr))
#   while X.filter(dWf)
#     dWr = step(mean, std, rng, T)
#     dWf = sum(real.(dWr))
#   end
#   X.Xpast .+= dWr

#   return dWf
# end

function (X::CARMAProcess)(t0::T, W0::T, Z0=nothing; method="CARMA", kwargs...) where {T<:Real}

  if method == "CARMA"
    if !isempty(X.carma.κ)
      # == Precalculate values == #
      sigma = X.arma.sigma .*  X.carma.α ./ sqrt.(-2.0 .*  X.carma.κ)

      # == Noise Process == #
      dist = CARMA_Process(X.carma.κ, X.carma.α, sigma, X.arma.mean, X.dt,zeros(ComplexF64,length(X.arma.phi)),length(X.arma.phi),X.Htime)
      return NoiseProcess{false}(float(t0),float(W0),Z0,dist,nothing; kwargs...)
    else
      return RealWienerProcess(float(t0),float(W0), Z0)
    end

  elseif method == "ARMA"
    return X.ARMAProcess(float(t0),float(W0), Z0)
  end

end