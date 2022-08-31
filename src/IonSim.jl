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

include("helpers.jl")
include("constants.jl")
include("ion/_include_ion.jl")
include("lasers.jl")
include("vibrational/_include_vibrational.jl")
include("ion_configurations.jl")
include("trap/_include_trap.jl")
include("ion_operators.jl")
include("hamiltonian/_include_hamiltonian.jl")
include("time_evolution.jl")

module analytical
include("analytic_functions.jl")
end

end  # main module
