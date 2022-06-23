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

include("constants/_include_constants.jl")

# First part of include ion (this means that ion/_include_ion.jl is not used)
include("ion/ion_properties.jl")
include("ion/ion.jl")
include("ion/ion_instance.jl")
include("ion/species/_include_species.jl")

include("lasers/_include_lasers.jl")
include("vibrational/_include_vibrational.jl")
include("ion_configurations/_include_ion_configurations.jl")
include("trap/_include_trap.jl")

# Second part of include ion (this means that ion/_include_ion.jl is not used)
include("ion/levels.jl")
include("ion/stark_shift.jl")
include("ion/zeeman_shift.jl")
include("ion/transitions.jl")

include("operators/_include_operators.jl")
include("hamiltonian/_include_hamiltonian.jl")
include("mathematical_functions/_include_mathematical_functions.jl")
include("convenience_functions.jl")
include("time_evolution.jl")

end  # main module
