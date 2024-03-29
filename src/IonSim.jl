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
include("vibrationalmodes.jl")
include("lasers.jl")
include("iontraps.jl")
include("chambers.jl")
include("operators.jl")
include("hamiltonians.jl")
include("timeevolution.jl")
include("species/_include_species.jl")

module analytical
include("analyticfunctions.jl")
end

end  # main module
