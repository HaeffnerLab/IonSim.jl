import StochasticDiffEq: calculate_step!, accept_step!

accept_step!(W::NoiseProcess,dt) = accept_step!(W, float(dt), nothing, nothing)
calculate_step!(W::NoiseProcess,dt) = calculate_step!(W, float(dt), nothing, nothing)

# == Calculate Noise Trajectory == #
function calculate_noise!(W::NoiseProcess, Nsteps::Integer)
    calculate_step!(W, W.dt)
    
    for step in 1:(length(W.u) == 1 ? (Nsteps-1) : Nsteps)
        accept_step!(W, W.dt)
    end

    return nothing
end