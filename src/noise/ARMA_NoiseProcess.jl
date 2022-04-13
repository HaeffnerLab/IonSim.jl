using DiffEqNoiseProcess
using LinearAlgebra: dot
using Random

@inline wiener_randn(rng::AbstractRNG,::Type{T}) where T = randn(rng,T)

# == Generators of Discrete AR/MA time series == #
struct AR_Process{T1,T2,T3,T4,T5,T6}
  phi::T1
  p::T2
  sigma::T3
  mean::T4
  c::T5
  filter::T6
end

mutable struct MA_Process{T1,T2,T3,T4,T5,T6}
  psi::T1
  q::T2
  sigma::T3
  mean::T4
  et::T5
  filter::T6
end

struct sinc_interp{T1,T2,T3,T4,T5}
  noise::T1
  dt::T2
  T::T3
  n::T4
  ω::T5
end

@inline function (X::AR_Process)(dW,W,dt,u,p,t,rng)
  T = (typeof(dW) <: AbstractArray) ? dW : typeof(dW)
  et = () -> X.sigma * wiener_randn(rng, T);

  # AR(p) Processes Step
  n = length(W.u)
  test = X.p >= n
  past = @fastmath X.c + dot(X.phi[1:(test ? n : end)], W.u[end:-1:end - (test ? n : X.p) + 1])

  new_val = past + et()
  while X.filter(new_val - X.mean) #Filter Function
    new_val = past + et()
  end

  return new_val - W.u[end]
end

@inline function (X::MA_Process)(dW,W,dt,u,p,t,rng)
  T = (typeof(dW) <: AbstractArray) ? dW : typeof(dW)
  et = () -> X.sigma * wiener_randn(rng, T)

  X.et[2:end] = X.et[1:end-1]
  X.et[1] = et()

  past = @fastmath X.mean + dot(X.psi[2:end],X.et[2:end])
  new_val = past + X.psi[1] * X.et[1]
  while X.filter(new_val - X.mean) #Filter Function
    X.et[1] = et()
    new_val = past + X.psi[1] * X.et[1]
  end

  return new_val - W.u[end]
end

@inline function WHITE_NOISE_BRIDGE(dW,arma,W,W0,Wh,q,h,u,p,t,rng)
  # Following correction uses Brownian Bridge
  # between W.dt values rather than nearest-neighbor
  Wcorrection = (1 - q) * W0
  h = W.dt
  q = (t - h*floor(t/h)) / h
  Wmean = q * W(h * ceil(t/h))[1] + (1 - q) * W(h * floor(t/h))[1]
  
  step = () -> @fastmath sqrt((1-q)*q*abs(h)) * arma.sigma * wiener_randn(rng, (typeof(dW) <: AbstractArray) ? dW : typeof(dW))

  dX = Wmean + step()
  while arma.filter(dX - arma.mean) #Filter Function
    dX = Wmean + step()
  end

  return dX - Wcorrection
end

@inline fourier_sinc(t, X) =  isapprox(mod(t, (X.T*X.dt)) + 1.0, 1.0) ? X.T : sin(X.n * X.ω * t) / sin(X.ω * t / 2.0)

@inline function (X::sinc_interp)(dW,W,dt,u,p,t,rng) # t::time at new step
 
  # == Sinc Function Interpolation of Noise == #
  tk = t .- Array(0:1:(X.T-1)).*X.dt
  tk[tk .< 0] .+= X.T*X.dt # periodic ϕ(-k) = ϕ(T-k)
  
  return @fastmath (1/X.T) * dot(X.noise, map(t->fourier_sinc(t, X), tk)) - W[end]
end

@inline function sinc_bridge(X,dW,W,W0,Wh,q,h,u,p,t,rng)
  
  # Brownian bridge correction term (added during interpolation)
  # iscontinuous(W) = true <- should fix to false
  Wc = (1 - q) * W0

  # == Sinc Function Interpolation of Noise == #
  tk = t .- Array(0:1:(X.T-1)).*X.dt
  tk[tk .< 0] .+= X.T*X.dt # periodic ϕ(-k) = ϕ(T-k)

  return @fastmath (1/X.T) * dot(X.noise, map(t->fourier_sinc(t, X), tk))  - Wc
end

# == ARMAProcess: Generates ARMA Coefficients and Sinc-interpolated NoiseProcess == #

mutable struct ARMAProcess
  psd::Function
  dt::Real
  Hfreq::Function
  Htime::Function
  normalize_sigma::Bool
  arma::ARMA
  window_length::Int
  timescale::Real

  function ARMAProcess(psd::Function, dt::Real; mean::Real=0.0, normalize_sigma::Bool=true, Hfreq::Union{Function,Nothing}=nothing, Htime::Union{Function,Nothing}=nothing, 
    window_length::Integer=2^10, coeff_length::Union{Integer,Nothing}=2^6, impulse_length::Union{Integer,Nothing}=nothing, timescale::Real=1.0,
    kwargs...)
    
    # == Temporal Filter == #
    #? Gaussian vs Uniform?
    if (Htime === nothing) 
      NoiseFilter = (x) -> false
    else
      if !(typeof(Htime(0.0)) <: Bool)
        throw(AssertionError("Htime(x) output must be of type Bool."))
      end
      NoiseFilter = (x) -> Htime(x)
    end

    # == Bandwidth Filter == #
    #? add DSP.jl compatability
    if !(Hfreq === nothing) 
      if !(typeof(Hfreq(0.0)) <: Real)
        throw(AssertionError("Hfreq(x) output must be of type Real."))
      end
      BandwidthFilter = (x) -> Hfreq(x)
      powerspectrum = (f) -> psd(f) * Hfreq(f)
    else
      BandwidthFilter = (x) -> 1.0
      powerspectrum = (f) -> psd(f)
    end
    
    # == Create ARMA Process == #
    dt = abs(dt)
    freq = fftfreq(window_length, 1 / (dt*timescale))
    ps = powerspectrum.(freq)

    # if (ps[1] == Inf) || (ps[1] == NaN) # est value from discrete fwd derivative
    #   Δ = freq[2] - freq[1]
    #   deriv = (-3 * ps[2] + 4 * ps[3] - ps[4])/(2*Δ)
    #   ps[1] = ps[2] + deriv * (freq[1] - freq[2])
    # end

    arma = ARMA(ps, freq; mean, normalize_sigma, coeff_length, impulse_length, dt, timescale)
    
    # # == Noise variance correction == #
    # n = length(arma.psd)÷2+1
    # psd_arma = map(f -> arma_psd(arma, f, dt * timescale), arma.freq[1:n])
    # arma.sigma = sqrt(sum((psd_arma[1:n] .- arma.psd[1:n]).^2)/n) / (sum(arma.psd[1:n])/n)

    # == Return ARMAProcess structure == #
    new(powerspectrum, float(dt), BandwidthFilter, NoiseFilter, normalize_sigma, arma, window_length, timescale)
  end
end

# == Sinc-interpolated NoiseProcess == #
function (X::ARMAProcess)(t0::T, W0::Union{T,Nothing}, Z0=nothing; method="AR", kwargs...) where {T<:Real}
  # == Useful == #
  arma = X.arma

  # psd_arma = arma_psd(arma, arma.freq[1], X.dt * X.timescale)
  # corr = abs(1.0 - (psd_arma[1] / arma.psd[1]))
  # sigma = arma.sigma * sqrt(1.0 + corr)
  # n = length(arma.psd)÷2+1
  # psd_arma = map(f -> arma_psd(arma, f, X.dt * X.timescale), arma.freq[1:n])
  # corr = sqrt(sum((1.0 .- psd_arma[1:n] ./ arma.psd[1:n]).^2)/n) # RMSPE of Power Spectrum
  # corr = sqrt(1.0 + corr)

  # == Initialize Random Start Value == #
  if W0 === nothing 
    # ensures random step on each redefinition
    step = () -> arma.sigma * randn()
    W0 = step()
    while X.Htime(W0)
      W0 = step()
    end
    W0 += arma.mean
  end

  # == Noise Process == #
  if method == "AR"
    dist = AR_Process(arma.phi,length(arma.phi), arma.sigma, arma.mean,arma.mean*(1-sum(arma.phi)),X.Htime)
  elseif (method == "Wold") || (method == "MA")
    dist = MA_Process(arma.psi,length(arma.psi), arma.sigma, arma.mean,[],X.Htime)
    dist.et = dist.sigma * randn(dist.q)
    dist.et[1] = (W0 - arma.mean)
  else
    err_msg = "Method must be either 'AR' (finite impulse response) or 'MA'/'Wold' (infinite impulse response)"
    @assert true err_msg
  end

  DiscreteProcess = NoiseProcess{false}(float(t0),float(W0),Z0,dist,(dW,W,W0,Wh,q,h,u,p,t,rng)->WHITE_NOISE_BRIDGE(dW,dist,W,W0,Wh,q,h,u,p,t,rng); kwargs...)
  DiscreteProcess.dt = X.dt
  return DiscreteProcess

  # ts_length = length(t0:X.dt:tend)
  # calculate_noise!(DiscreteProcess, ts_length < 1024 ? 1024 : ts_length) #? Consider inplace alternative
  # N = length(DiscreteProcess.u)
  # dist = sinc_interp(DiscreteProcess.u, DiscreteProcess.dt, N, (N-1.0)/2.0, 2pi/(N*X.dt))

  # return NoiseProcess{false}(float(t0),float(W0),Z0,dist,(dW,W,W0,Wh,q,h,u,p,t,rng) -> sinc_bridge(dist,dW,W,W0,Wh,q,h,u,p,t,rng); kwargs...)
end