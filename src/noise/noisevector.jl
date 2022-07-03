import DiffEqNoiseProcess: AbstractNoiseProcess
import DiffEqNoiseProcess: interpolate!

mutable struct NoiseVector{T, N, Tt, T2, T3, ZType, inplace} <:
               AbstractNoiseProcess{T, N, Vector{T2}, inplace}
    t::Vector{Tt}
    u::Vector{T2}
    W::Vector{T2}
    Z::ZType
    curt::Tt
    curW::T2
    curZ::T3
    dt::Tt
    dW::T2
    dZ::T3
    step_setup::Bool
    reset::Bool
    linear::Bool
end

function NoiseVector(t, W, Z = nothing; reset = true, linear = false)
    val = W[1]
    curt = t[1]
    dt = t[1]
    curW = copy(val)
    dW = copy(val)
    if Z == nothing
        curZ = nothing
        dZ = nothing
    else
        curZ = copy(Z[1])
        dZ = copy(Z[1])
    end
    (typeof(val) <: AbstractArray && !(typeof(val) <: SArray)) ? iip = true : iip = false
    return NoiseVector{
        typeof(val),
        ndims(val),
        typeof(dt),
        typeof(dW),
        typeof(dZ),
        typeof(Z),
        iip
    }(
        t,
        W,
        W,
        Z,
        curt,
        curW,
        curZ,
        dt,
        dW,
        dZ,
        true,
        reset,
        linear
    )
end

(W::NoiseVector)(t) = interpolate!(W, t)
(W::NoiseVector)(u, p, t) = interpolate!(W, t)
(W::NoiseVector)(out1, out2, u, p, t) = interpolate!(out1, out2, W, t)

# == Interface Functions == #
# https://github.com/SciML/DiffEqNoiseProcess.jl/blob/c48cdce099cece1edbd8f99da960bc67e3c2c4ca/src/noise_interfaces/noise_grid_interface.jl

function interpolate!(W::NoiseVector, t)
    ts, timeseries, timeseries2 = W.t, W.W, W.Z
    sign(W.dt) * t > sign(W.dt) * ts[end] && error(
        "Solution interpolation cannot extrapolate past the final timepoint. Build a longer NoiseVector to cover the integration."
    )
    sign(W.dt) * t < sign(W.dt) * ts[1] && error(
        "Solution interpolation cannot extrapolate before the first timepoint. Build a longer NoiseVector to cover the integration."
    )
    tdir = sign(ts[end] - ts[1])

    if t isa Union{Rational, Integer}
        @inbounds i = searchsortedfirst(ts, t, rev = tdir < 0) # It's in the interval ts[i-1] to ts[i]
    else
        @inbounds i = searchsortedfirst(ts, t - tdir * 10eps(typeof(t)), rev = tdir < 0)
    end

    @inbounds if (t isa Union{Rational, Integer} && ts[i] == t) ||
                 (isapprox(t, ts[i]; atol = 100eps(typeof(t)), rtol = 100eps(t)))
        val1 = timeseries[i]
        timeseries2 !== nothing ? val2 = timeseries2[i] : val2 = nothing
    elseif ts[i - 1] == t # Can happen if it's the first value!
        val1 = timeseries[i - 1]
        timeseries2 !== nothing ? val2 = timeseries2[i - 1] : val2 = nothing
    else
        if W.linear
            dt = ts[i] - ts[i - 1]
            Θ = (t - ts[i - 1]) / dt
            val1 = linear_interpolant(Θ, dt, timeseries[i - 1], timeseries[i])
            timeseries2 !== nothing ?
            val2 = linear_interpolant(Θ, dt, timeseries2[i - 1], timeseries2[i]) :
            val2 = nothing
        else
            val1 = timeseries[i - 1]
            timeseries2 !== nothing ? val2 = timeseries2[i - 1] : val2 = nothing
        end
    end
    return isa(vals2, Nothing) ? val1 : (val1, val2)
end

function interpolate!(out1, out2, W::NoiseVector, t)
    ts, timeseries, timeseries2 = W.t, W.W, W.Z
    sign(W.dt) * t > sign(W.dt) * (ts[end] + 10 * sign(W.dt) * eps(typeof(t))) && error(
        "Solution interpolation cannot extrapolate past the final timepoint. Build a longer NoiseVector to cover the integration."
    )
    sign(W.dt) * t < sign(W.dt) * (ts[1] - 10 * sign(W.dt) * eps(typeof(t))) && error(
        "Solution interpolation cannot extrapolate before the first timepoint. Build a longer NoiseVector to cover the integration."
    )
    tdir = sign(ts[end] - ts[1])

    if t isa Union{Rational, Integer}
        @inbounds i = searchsortedfirst(ts, t, rev = tdir < 0) # It's in the interval ts[i-1] to ts[i]
    else
        @inbounds i = searchsortedfirst(ts, t - tdir * 10eps(typeof(t)), rev = tdir < 0)
    end

    @inbounds if (t isa Union{Rational, Integer} && ts[i] == t) ||
                 (isapprox(t, ts[i]; atol = 100eps(typeof(t)), rtol = 100eps(t)))
        copyto!(out1, timeseries[i])
        timeseries2 !== nothing && copyto!(out2, timeseries2[i])
    elseif ts[i - 1] == t # Can happen if it's the first value!
        copyto!(out1, timeseries[i - 1])
        timeseries2 !== nothing && copyto!(out2, timeseries2[i - 1])
    else
        if W.linear
            dt = ts[i] - ts[i - 1]
            Θ = (t - ts[i - 1]) / dt
            linear_interpolant!(out1, Θ, dt, timeseries[i - 1], timeseries[i])
            timeseries2 !== nothing &&
                linear_interpolant!(out2, Θ, dt, timeseries2[i - 1], timeseries2[i])
        else
            out1 .= timeseries[i - 1]
            timeseries2 !== nothing && out2 .= timeseries2[i - 1]
        end
    end
end
