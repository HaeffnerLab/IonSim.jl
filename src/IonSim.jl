module IonSim

using QuantumOptics
import OrderedCollections: OrderedDict

export QuantumOptics, OrderedDict
export analytical

Base.copy(x::T) where T = T([getfield(x, k) for k ∈ fieldnames(T)]...)

include("constants.jl")
using .PhysicalConstants
include("ions.jl")
include("vibrational_modes.jl")
include("lasers.jl")
include("ion_configurations.jl")
include("traps.jl")
include("operators.jl")                           
include("hamiltonians.jl")
include("time_evolution.jl")

module analytical
    include("analytic_functions.jl")
end

end  # main module
