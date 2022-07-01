
# == ARMAProcess: Generates ARMA Coefficients == #
"""
bar(x[, y])

Compute the Bar index between `x` and `y`.

If `y` is unspecified, compute the Bar index between all pairs of columns of `x`.

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
mutable struct ARMAProcess
	arma::ARMA
	powerspectrum::Union{Function,Vector}
	Hfreq::Function
	Htime::Function
	τ::Real
	timescale::Real

	@inline function ARMAProcess(psd::Function, dt::Real, timescale::Real=1.0; mean::Real=0.0, normalize_sigma::Bool=true, Hfreq::Union{Function,Nothing}=nothing, Htime::Union{Function,Nothing}=nothing,
	  window_length::Integer=2^10, coeff_length::Union{Integer,Nothing}=2^6, impulse_length::Union{Integer,Nothing}=nothing, rtol::Real=1e-6)

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
		BandwidthFilter = Hfreq
		powerspectrum = (f) -> psd(f) * Hfreq(f)
	  else
		BandwidthFilter = (x) -> 1.0
		powerspectrum = (f) -> psd(f)
	  end

	  # == Create ARMA Process == #
	  dt = abs(dt)
	  freq = fftfreq(2window_length, 1 / (dt*timescale))
	  ps = powerspectrum.(abs.(freq))
	  arma = ARMA(ps, freq, float(dt), timescale; mean, normalize_sigma, coeff_length, impulse_length, rtol)

	  # == Return ARMAProcess structure == #
	  new(arma, powerspectrum, BandwidthFilter, NoiseFilter, float(dt), timescale)
	end
end

@inline function (X::ARMAProcess)(t0::T, W0::T=0.0, Z0=nothing; response="AR", kwargs...) where {T<:Real}
	# == Useful == #
	arma = X.arma

	# == Noise Process == #
	if response == "AR"
	  dist = AR_Process(arma.phi,length(arma.phi), arma.sigma, arma.mean, arma.mean*(1-sum(arma.phi)),X.arma.dt,float(t0) + X.arma.dt,X.Htime, zeros(length(arma.phi)))
	  for i in 1:min(256,2*dist.p)
		mean = dist.c + dot(dist.phi, dist.past)
		pushfirst!(dist.past, 0.0); pop!(dist.past)
		dist.past[1] = mean + dist.sigma * randn()
	  end
	  dist.past[1] = W0
	  bridge = AR_STEP_BRIDGE
	elseif (response == "Wold") || (response == "MA")
	  dist = MA_Process(arma.psi,length(arma.psi), arma.sigma, X.arma.dt, float(t0) + X.arma.dt, arma.mean,zero(arma.psi),X.Htime)
	  dist.et .= dist.sigma * randn(dist.q)
	  dist.et[1] = (W0 - arma.mean)
	  bridge = MA_STEP_BRIDGE
	else
	  err_msg = "response must be either 'AR' (finite impulse response) or 'MA'/'Wold' (infinite impulse response)"
	  @assert true err_msg
	end

	return NoiseProcess{false}(float(t0),float(W0),Z0,dist,(dW,W,W0,Wh,q,h,u,p,t,rng)->bridge(dist,dW,W,W0,Wh,q,h,u,p,t,rng); rswm=RSWM(adaptivealg=:RSwM3), kwargs...)
end

# == CARMAProcess struct & function == #
"""
bar(x[, y])

Compute the Bar index between `x` and `y`.

If `y` is unspecified, compute the Bar index between all pairs of columns of `x`.

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
struct CARMAProcess{T1,T2,T3,T4,ARMA,T5}
	αr::T1
	λr::T1
	σr::T1
	psd::T2
	acf::T3
	mean::T4
	arma::ARMA
	Htime::T5
	timescale::T4

	@inline function CARMAProcess(arma::ARMA; Htime::Function = x -> false, fit_args::Union{Dict,Nothing}=nothing)

		# == Check for non-zero AR coefficients == #
		if any(abs.(arma.phi) .> eps())
			# == continous model poles == #
			carma_λ = ARMAtoCARMAPoles(arma.ar_roots, arma.dt)

			# == Fit Discrete to Continuous ACF Model == #
			#! TODO: add warn statement if poor fit; use statistical measurement for goodness of fit
			carma_α, carma_sigma, carma_psd, carma_acf = ARMAtoCARMASolver(arma.dt, arma.acf, arma.phi, arma.psi, carma_λ, max(arma.sigma, arma.ma_sigma); fit_args=fit_args)

			# == Precalculate part of carma sigma == #
			carma_sigma = @. carma_sigma * carma_α / sqrt(-2.0 *  carma_λ)

		else # zero coefficients means white noise model
			carma_λ = carma_α = carma_sigma = carma_acf = carma_psd = []
		end

		# == Return CARMAProcess structure == #
		new{typeof(carma_α),typeof(carma_psd), typeof(carma_acf),typeof(arma.mean),typeof(arma),typeof(Htime)}(carma_α, carma_λ, carma_sigma, carma_psd, carma_acf, arma.mean, arma, Htime, arma.timescale)
	end
end

# may take ARMAProcess as input
CARMAProcess(arma::ARMAProcess; kwargs...) = CARMAProcess(arma.arma; Htime=arma.Htime, kwargs...)

@inline function (X::CARMAProcess)(t0::T, W0::T, Z0=nothing; kwargs...) where {T<:Real}
	# past_values = Dict(float(t0) => tmp) # Dict to cache values
	dist = CARMA_Process(X.αr, X.λr, X.mean, X.σr, float(t0), X.Htime)
	# initialize noise source memory
	# for i in 1:min(256,length(X.λr))
	  # carma_meanstd!(dist.mean,dist.std,dist.Xpast,dist.sigma,dist.λ,X.arma.dt)
	  # carma_step!(dist.Xpast,dist.mean,dist.std)
	# end
	return NoiseProcess{false}(float(t0),float(W0),Z0,dist,(dW,W,W0,Wh,q,h,u,p,t,rng) -> CARMABridge(dist,dW,W,W0,Wh,q,h,u,p,t,rng); rswm=RSWM(adaptivealg=:RSwM2), kwargs...)
end

# == GaussianProcess struct & function == #
"""
bar(x[, y])

Compute the Bar index between `x` and `y`.

If `y` is unspecified, compute the Bar index between all pairs of columns of `x`.

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
mutable struct GaussianProcess
	ARMAProcess::ARMAProcess
	CARMAProcess::Any
	continuous::Bool
	dt::Real
	timescale::Real

	@inline function GaussianProcess(psd::Function, dt::Real, timescale::Real=1.0; mean::Real=0.0, normalize_sigma::Bool=true, Hfreq::Union{Function,Nothing}=nothing, Htime::Union{Function,Nothing}=nothing,
		window_length::Integer=2^10, coeff_length::Union{Integer,Nothing}=2^7, impulse_length::Union{Integer,Nothing}=nothing, rtol::Real=1e-6, continuous::Bool=false, fit_args::Union{Dict,Nothing}=nothing)

		# == Discrete ARMA Model coefficients == #
		arma = ARMAProcess(psd, dt, timescale; mean, normalize_sigma, Hfreq, Htime,
		window_length, coeff_length, impulse_length, rtol=rtol)

		# == Continuous ARMA Model coefficients == #
		# avoid initializing if not needed
		if continuous
		  carma = CARMAProcess(arma; fit_args=fit_args)
		  if isempty(carma.αr)
			  carma = WienerProcess
		  end
		else
		  carma = nothing
		end

		# == Return ARMAProcess structure == #
		new(arma,carma,continuous,dt,timescale)
	end
end

@inline function (X::GaussianProcess)(t0::T, W0::T, Z0=nothing; continuous::Bool=false, kwargs...) where {T<:Real}
	if continuous
		try
		  return X.CARMAProcess(float(t0),float(W0), Z0; kwargs...)
		catch
		  X.continuous = true
		  return X.CARMAProcess(float(t0),float(W0), Z0; kwargs...)
		end
	else #discrete model
		return X.ARMAProcess(float(t0),float(W0), Z0; kwargs...)
	end
end

# == Setting Property on GaussianProcess == #
function Base.setproperty!(X::GaussianProcess, s::Symbol, v::Tv) where {Tv}
  if (s == :continuous)
	if v && isa(X.CARMAProcess,Nothing)
	  vnew = CARMAProcess(X.ARMAProcess)
	  if isempty(vnew.αr)
		  vnew = WienerProcess
	  end
	elseif !v
	  vnew = nothing
	end
	Core.setproperty!(X, :CARMAProcess, vnew)
	Core.setproperty!(X, s, v)
	return
  end

  # Other parameters involve recalculating ARMA
  discrete = X.ARMAProcess
  arma = X.ARMAProcess.arma

  if (s == :dt)
	vnew = ARMAProcess(discrete.powerspectrum, s, arma.timescale; arma.mean, arma.normalize_sigma, discrete.Hfreq, discrete.Htime,
		arma.window_length, arma.coeff_length, arma.impulse_length, arma.arma_rtol)
  elseif (s == :timescale)
	vnew = ARMAProcess(discrete.powerspectrum, arma.dt, s; arma.mean, arma.normalize_sigma, discrete.Hfreq, discrete.Htime,
	arma.window_length, arma.coeff_length, arma.impulse_length, arma.arma_rtol)
  else
	error("Property cannot be changed. Recalculate GaussianProcess instead.")
  end

  Core.setproperty!(X, s, v)
  Core.setproperty!(X, :ARMAProcess, vnew)
  return
end

#! TODO: Create ARMAtoCARMA function that takes a pre-made ARMAProcess