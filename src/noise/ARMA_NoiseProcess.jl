using DiffEqNoiseProcess
using LinearAlgebra: dot
using Random
using PyCall

@inline wiener_randn(rng::AbstractRNG, ::Type{T}) where {T} = randn(rng, T)
function arma_step!(x, mean, sigma, rng, T)
    return x[1] = mean + sigma * wiener_randn(rng, T)
end

# == Generators of Discrete AR/MA time series == #
mutable struct AR_Process{T1, T2, T3, T4, T5}
    phi::T1
    p::T2
    sigma::T3
    mean::T3
    c::T3
    dt::T3
    t::T3
    filter::T4
    past::T5
end

mutable struct MA_Process{T1, T2, T3, T4, T5}
    psi::T1
    q::T2
    sigma::T3
    mean::T3
    dt::T3
    t::T3
    et::T4
    filter::T5
end

@inline function (X::AR_Process)(dW, W, dt, u, p, t, rng)
    T = (typeof(dW) <: AbstractArray) ? dW : typeof(dW)

    # == derivative of step function == #
    ϵ = t isa Union{Rational, Integer} ? 0 : 100eps(typeof(t))
    t += dt
    if t < X.t + ϵ
        return 0 #dW = 0 between steps
    end

    # == AR(p) Step == #
    while t >= X.t + ϵ
        mean = X.c + dot(X.phi, X.past)
        pushfirst!(X.past, 0.0)
        pop!(X.past)

        arma_step!(X.past, mean, X.sigma, rng, T)
        while X.filter(X.past[1] - X.mean) #Filter Function
            arma_step!(X.past, mean, X.sigma, rng, T)
        end

        X.t += X.dt
    end

    return X.past[1] * max(dt, X.dt)
end

@inline function (X::MA_Process)(dW, W, dt, u, p, t, rng)
    T = (typeof(dW) <: AbstractArray) ? dW : typeof(dW)

    # == derivative of step function == #
    ϵ = t isa Union{Rational, Integer} ? 0 : 100eps(typeof(t))
    t += dt
    if t < X.t + ϵ
        return 0 #dW = 0 between steps
    end

    # == MA(∞) Step == #
    Xt = 0
    while t >= X.t + ϵ
        pushfirst!(X.et, 0.0) #add new first value
        pop!(X.et) #remove last value
        arma_step!(X.et, 0.0, X.sigma, rng, T)
        Xt = X.mean + dot(X.psi, X.et)

        while X.filter(Xt - X.mean) #Filter Function
            arma_step!(X.et, 0, X.sigma, rng, T)
            Xt = X.mean + dot(X.psi, X.et)
        end

        X.t += X.dt
    end
    return Xt * max(dt, X.dt)
end

# == Bridge Function == #
@inline function AR_STEP_BRIDGE(X, dW, W, W0, Wh, q, h, u, p, t, rng)
    # setup_next_step! correction
    if isapprox(W0, W.curW) # line 97
        return AR_STEP_BRIDGE(X, dW, W, 1, Wh, q, h, u, p, t, rng)
    end

    if t isa Union{Rational, Integer}
        i = searchsortedlast(W.t, h * floor(t / h))
    else
        i = searchsortedlast(W.t, h * floor(t / h) + eps(typeof(t)))
    end
    return W[i] - (1 - q) * W0
end

# reject_step! correction
# https://github.com/SciML/DiffEqNoiseProcess.jl/blob/master/src/noise_interfaces/noise_process_interface.jl
@inline function AR_STEP_BRIDGE(X, dW, W, W0::Int, Wh, q, h, u, p, t, rng)

    # line 203 correction
    # q = dtnew / W.dt -> h = W.dt (not dtnew)
    if isempty(W.S₂) && isapprox(h, q * W.dt) # if h=dtnew
        h = W.dt
    end

    tnew = t + q * h
    tfloor = (
        X.dt * floor(round(t / X.dt, digits = 8)),
        X.dt * floor(round(tnew / X.dt, digits = 8))
    )
    ϵ = t isa Union{Rational, Integer} ? 0 : sqrt(eps(typeof(t)))

    while tnew < (X.t - X.dt) + ϵ
        X.t -= X.dt
        popfirst!(X.past)
        push!(X.past, 0.0)
        # X.past .= [X.past[2:end]; 0.0]
    end

    # both in same step
    if isapprox(tfloor[1], tfloor[2])
        dWnew = 0.0
        # both in different steps
    else
        dWnew = X.past[1] * min(q * h, X.dt)
    end

    return dWnew
end

@inline function MA_STEP_BRIDGE(X, dW, W, W0, Wh, q, h, u, p, t, rng)
    # setup_next_step! correction
    if isapprox(W0, W.curW) # line 97
        return MA_STEP_BRIDGE(X, dW, W, 1, Wh, q, h, u, p, t, rng)
    end

    if t isa Union{Rational, Integer}
        i = searchsortedlast(W.t, h * floor(t / h))
    else
        i = searchsortedlast(W.t, h * floor(t / h) + eps(typeof(t)))
    end
    return W[i] - (1 - q) * W0
end

# reject_step! correction
# https://github.com/SciML/DiffEqNoiseProcess.jl/blob/master/src/noise_interfaces/noise_process_interface.jl
@inline function MA_STEP_BRIDGE(X, dW, W, W0::Int, Wh, q, h, u, p, t, rng)

    # line 203 correction
    # q = dtnew / W.dt -> h = W.dt (not dtnew)
    if isempty(W.S₂) && isapprox(h, q * W.dt) # if h=dtnew
        h = W.dt
    end

    tnew = t + q * h
    tfloor = (
        X.dt * floor(round(t / X.dt, digits = 8)),
        X.dt * floor(round(tnew / X.dt, digits = 8))
    )
    ϵ = t isa Union{Rational, Integer} ? 0 : sqrt(eps(typeof(t)))

    while tnew < (X.t - X.dt) + ϵ
        X.t -= X.dt
        popfirst!(X.et) # remove last innovation
        push!(X.et, 0.0) # add zero to end
    end

    # both in same step
    if isapprox(tfloor[1], tfloor[2])
        dWnew = 0.0
        # both in different steps
    else
        dWnew = (X.mean + dot(X.psi, X.et)) * min(q * h, X.dt)
    end

    return dWnew
end

# @inline function WHITE_NOISE_BRIDGE(dW,arma,W,W0,Wh,q,h,u,p,t,rng)
#   # Following correction uses Brownian Bridge
#   # between W.dt values rather than nearest-neighbor
#   # Wcorrection = (1 - q) * W0
#   # h = W.dt
#   # q = (t - h*floor(t/h)) / h
#   # Wmean = q * W(h * ceil(t/h))[1] + (1 - q) * W(h * floor(t/h))[1]
#   mean = q * Wh
#   sigma = arma.sigma * sqrt((1-q)*q*abs(h))

#   T = (typeof(dW) <: AbstractArray) ? dW : typeof(dW)
#   Xf = et(mean,sigma,rng,T)
#   while arma.filter(step - arma.mean) #Filter Function
#     Xf = et(mean,sigma,rng,T)
#   end

#   return Xf
# end

# == potentially useful sinc interpolation stuff == #
# @inline fourier_sinc(t, X) =  isapprox(mod(t, (X.T*X.dt)) + 1.0, 1.0) ? X.T : sin(X.n * X.ω * t) / sin(X.ω * t / 2.0)

# struct sinc_interp{T1,T2,T3,T4,T5}
#   noise::T1
#   dt::T2
#   T::T3
#   n::T4
#   ω::T5
# end

# @inline function (X::sinc_interp)(dW,W,dt,u,p,t,rng) # t::time at new step

#   # == Sinc Function Interpolation of Noise == #
#   tk = t .- Array(0:1:(X.T-1)).*X.dt
#   tk[tk .< 0] .+= X.T*X.dt # periodic ϕ(-k) = ϕ(T-k)

#   return @fastmath (1/X.T) * dot(X.noise, map(t->fourier_sinc(t, X), tk)) - W[end]
# end

# @inline function sinc_bridge(X,dW,W,W0,Wh,q,h,u,p,t,rng)

#   # Brownian bridge correction term (added during interpolation)
#   # iscontinuous(W) = true <- should fix to false
#   Wc = (1 - q) * W0

#   # == Sinc Function Interpolation of Noise == #
#   tk = t .- Array(0:1:(X.T-1)).*X.dt
#   tk[tk .< 0] .+= X.T*X.dt # periodic ϕ(-k) = ϕ(T-k)

#   return @fastmath (1/X.T) * dot(X.noise, map(t->fourier_sinc(t, X), tk))  - Wc
# end
