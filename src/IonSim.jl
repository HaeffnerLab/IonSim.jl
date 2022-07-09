module IonSim

using QuantumOptics
using QuantumOpticsBase: BASES_CHECK
import OrderedCollections: OrderedDict

#=
This is a band-aid for https://github.com/HaeffnerLab/IonSim.jl/issues/90
For now, we just turn of bases checks. This means that if users feed an initial state and 
Hamiltonian (and/or collapse operators) to the solvers that have incompatbile bases an opaque
exception will be thrown. This will be especially annoying when e.g. vibrational modes are
misordered in a tensor product, since the correct order is implicilty enforced atm by the 
way hamiltonian is constructed.
TODO: wrap QO solvers with an IonSim function that handles things like bases checks.
=#
BASES_CHECK.x = false

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

module analytical
include("analytic_functions.jl")
end

end  # main module
