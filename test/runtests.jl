# Might be a bad idea, but it's convenient! If there's a more julia-esque idiom
# for test fixtures, that's much preferred.
using IonSim.Properties: loadfromconfig
CA40_PROPERTIES = loadfromconfig("../configs/ions/ca40.yaml")

include("test_analytic_functions.jl")
include("test_constants.jl")
include("test_operators.jl")
include("test_lasers.jl")
include("test_vibrational_modes.jl")
include("test_ions.jl")
include("test_ion_traps.jl")
include("test_chambers.jl")
include("test_hamiltonians.jl")
include("test_dynamics.jl")
