module GaussianNoiseProcess

include("ARMA_Structure.jl")
include("ARMA_NoiseProcess.jl")
include("CARMA_NoiseProcess.jl")
include("ARMA_StepInterface.jl")

export ARMA

export GaussianProcess

export ARMAProcess

export CARMAProcess

export calculate_step!, accept_step!, calculate_noise!

end