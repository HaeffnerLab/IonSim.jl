module IonSim

using QuantumOptics
import OrderedCollections: OrderedDict

export OrderedDict, analytical
# Export some commonly used QuantumOptics.jl functions
export embed,
    ⊗,
    tensor,
    ⊕,
    dagger,
    normalize,
    normalize!,
    expect,
    tr,
    ptrace,
    tracenorm,
    tracedistance,
    entropy_vn,
    fidelity,
    diagonaljumps,
    dm,
    exp,
    norm

# used for copying composite types
Base.copy(x::T) where {T} = T([getfield(x, k) for k in fieldnames(T)]...)
"""
    IonSimBasis
An abstract type for specialized bases, which are unique to IonSim.
"""
abstract type IonSimBasis <: Basis end
export IonSimBasis

include("constants.jl")
include("ions.jl")
include("vibrational_modes.jl")
include("lasers.jl")
include("ion_configurations.jl")
include("traps.jl")
include("operators.jl")
include("hamiltonians.jl")
include("time_evolution.jl")
include("species/_include_species.jl")

module GaussianNoiseProcess

include("noise/ARMA_Structure.jl")
include("noise/ARMA_NoiseProcess.jl")
include("noise/CARMA_NoiseProcess.jl")
include("noise/ARMA_StepInterface.jl")
export ARMA, ARMAProcess, CARMAProcess
export calculate_step!, accept_step!, calculate_noise!

end

module analytical
include("analytic_functions.jl")
end

end  # main module
