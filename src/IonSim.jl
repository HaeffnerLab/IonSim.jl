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
include("noise/noisevector.jl")
include("lasers.jl")
include("ion_configurations.jl")
include("traps.jl")
include("operators.jl")
include("hamiltonians.jl")
include("time_evolution.jl")
include("species/_include_species.jl")


module noise
using FFTW
using LinearAlgebra: dot
using IterativeSolvers: lsqr, lsqr!
using Random
using StaticArrays
using PyCall
np = pyimport("numpy")
so = pyimport("scipy.optimize")

include("noise/ARMA_NoiseProcess.jl")
include("noise/ARMA_Structure.jl")
include("noise/CARMA_NoiseProcess.jl")
include("noise/CARMA_Structure.jl")
include("noise/Gaussian_types.jl")
include("noise/Gaussian_interface.jl")
export GaussianProcess, CARMAProcess, ARMAProcess, ARMA, calculate_noise!
# export calculate_step!, accept_step!, calculate_noise!
end

module analytical
include("analytic_functions.jl")
end

end  # main module
