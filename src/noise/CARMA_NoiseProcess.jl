# == Helper Functions == #

# @inline carma_randn(rng::AbstractRNG,::Type{T}) where T = randn(rng,T) + 1im * randn(rng,T)

@inline function carma_step!(x::T, mean::T, std::T, rng::AbstractRNG, Trng) where {T}
  rndm = randn(rng,Trng) #carma_randn(rng, Trng)
  @simd for i in eachindex(x)
    @inbounds x[i] = mean[i] + rndm * std[i]
  end
  return real(sum(x))
end

function carma_step!(x::T, mean::T, std::T) where {T}
  @simd for i in eachindex(x)
    @inbounds x[i] = mean[i] + randn() * std[i]
  end
end

@inline function step_n_filter!(x::T, filter::Function, mean::T, std::T, rng::AbstractRNG, Trng) where {T}
  Wf = carma_step!(x, mean, std, rng, Trng)
  while filter(Wf)
    Wf = carma_step!(x, mean, std, rng, Trng)
  end
  return Wf
end

@inline function carma_meanstd!(mean, std, x, sigma, λ, dt)
  @simd for i in eachindex(x)
    @inbounds mean[i] = x[i] * exp(λ[i]*dt)
    @inbounds std[i] = sigma[i] * sqrt(1.0 - exp(λ[i]*2.0dt)) # X.α[i] / sqrt(-2.0 * X.λ[i])
  end
end

@inline function carma_bridge_meanstd!(mean,std,Wi,Wh,sigma,λ,q,h)
  tmp = @. (sinh(λ * (1-q)*h), sinh(λ * q*h), sinh(λ * h))
  @simd for i in eachindex(mean)
      @inbounds mean[i] = (Wi[i] * tmp[1][i] + Wh[i] * tmp[2][i]) / tmp[3][i]
      @inbounds std[i] = sigma[i] * sqrt(-2.0 * tmp[2][i] * tmp[1][i] / tmp[3][i])
      # correct for redef of sigma in GaussianProcess
      # sigma = X.arma.sigma *  X.carma.α / sqrt(2.0 *  X.carma.λ)
  end
end

# == CARMA Stepping Process == #
mutable struct CARMA_Process{T1<:SArray,T2<:AbstractFloat,T3<:AbstractArray,T4<:Bool,T5<:Function}
  λ::T1
  sigma::T1
  μ::T2
  t::T2
  Xpast::T3
  mean::T3
  std::T3
  filter::T5

  function CARMA_Process(α, λ, μ, σ, t, filter)
    N = length(λ); Tel = eltype(λ)
    λ = SVector{N}(λ)
    sigma = SVector{N}(σ)
    tmp =zeros(Tel,N)
    return new{typeof(λ),typeof(μ),typeof(tmp),Bool,typeof(filter)}(λ,sigma,μ,t,tmp,zero(tmp),zero(tmp),filter)
  end
end

@inline function (X::CARMA_Process)(dW,W,dt,u,p,t,rng)
  T = (typeof(dW) <: AbstractArray) ? dW : typeof(dW) #? Still necessary?

  # == calculate step == #
  # https://julialang.org/blog/2013/09/fast-numeric/
  X.t = t+dt
  carma_meanstd!(X.mean,X.std,X.Xpast,X.sigma,X.λ,dt)
  return step_n_filter!(X.Xpast,X.filter,X.mean,X.std,rng,T) * dt

  # == past memory problem == #
  # X.Xpast[t+dt] = Wf[1] # past dictionary expensive (?)
end

# == CARMA Bridge == #
@inline function CARMABridge(X,dW,W,W0,Wh,q,h,u,p,t,rng)
  # see https://arxiv.org/pdf/1011.0548.pdf
  # pg 20 (space-time transform)
  # t ∈ [0, T]

  # == setup_next_carma_step!/reject_carma_step! correction == #
  if isapprox(W0,W.curW) || isa(W0,Int) #line 97 / 283
    # W0 = 0.0
    # solves dtmin crashing problem
    return 0.0
  else
    # == Bridging values == #
    t -= q * h
    q = q * h / (X.t - t)
    h = (X.t - t)
  end

  Wi = X.mean; fill!(Wi,0.0)
  Wf = X.Xpast
  # t = (t > eps()) ? searchkeys(reject_step ? t : t - q * h, X.Xpast) : 0.0
  # Wi = X.mean
  # Wi .= reject_step ? zeros(ComplexF64, length(X.λ)) : X.Xpast[t0]
  # Wh = X.std
  # Wh .= X.Xpast[searchkeys(t0 + h, X.Xpast)]  t0 = (t > eps()) ? (reject_step ? t : t - q * h) : 0.0

  # == Calculate Mean and STD of Bridge == #
  carma_bridge_meanstd!(X.mean,X.std,Wi,Wf,X.sigma,X.λ,q,h)

  # == Step Process == #
  T = (typeof(dW) <: AbstractArray) ? dW : typeof(dW) #? Still necessary?
  Wbridge = step_n_filter!(Wf,X.filter,X.mean,X.std,rng,T)

  # == add Bridge values to history == #
  # X.Xpast[t0 + (q * h)] = Wbridge[1]

  # == Bridge fix: correction to Brownian Bridge mean flow == #
  return (Wbridge * q * h) + q * W0 #- (1 - q) * W0 # line 445 (link above)
end

# # == reject_carma_step! fix == #
# #see https://github.com/SciML/DiffEqNoiseProcess.jl/blob/master/src/noise_interfaces/noise_process_interface.jl
# @inline function CARMABridge(X,dW,W,W0::Int,Wh,q,h,u,p,t,rng)
#   T = (typeof(dW) <: AbstractArray) ? dW : typeof(dW) #? Still necessary?

#   # line 203 correction
#   # q = dtnew / W.dt -> h = W.dt (not dtnew)
#   if isempty(W.S₂) && isapprox(h, q*W.dt) # if h=dtnew
#     h = W.dt
#   end

#   # == Bridging values == #
#   q = q * h / (X.t - t)
#   h = (X.t - t)
#   Wi = X.mean; fill!(Wi,0.0)
#   Wh = X.Xpast

#   # == Calculate Mean and STD of Bridge == #
#   carma_bridge_meanstd!(X.mean,X.std,Wi,Wh,X.sigma,X.λ,q,h)

#   # == Step Process == #
#   X.t = t + q*h
#   Wbridge = step_n_filter!(Wh,X.filter,X.mean,X.std,rng,T)

#   return  Wbridge * q * h
# end

# @inline function searchkeys(t, dict)
#   if haskey(dict, t)
#     return t
#   else
#     klist = collect(keys(dict))
#     return klist[findmin(abs, klist .- t)[2]]
#     # key = sort(collect(keys(dict)))
#     # index = searchsortedlast(key, t + (t isa Union{Rational,Integer} ? 0 : atol))
#     # return key[index]
#   end
# end