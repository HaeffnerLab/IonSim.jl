
# == Calculate Noise Trajectory == #
function calculate_noise!(W::NoiseProcess, tspan)
    dtlist = diff(tspan)
    if any(dtlist .!= dtlist[1])
        map((t,dt) -> (W.dt=dt; W(t)), tspan, [dtlist; dtlist[end]])
    else
        W.dt = dtlist[1]
        map(t -> W(t), tspan)
    end

    return (W.t, [0.0; diff(W.u)./diff(W.t)])
end