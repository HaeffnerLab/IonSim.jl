# This file contains functionality for tracking and manipulating ion properties.
module Properties  # Open to better names for this module.

import YAML
using LinearAlgebra: dot
using OrderedCollections

export IonProperties, loadfromconfig

function roundnearesthalf(value::Number)::Rational
    halves::Integer = 2 * value
    return halves // 2
end


"""
IonProperties type.
**Required keywords**
* `mass`: Mass of ion in kg
* `charge`: Charge of ion in units of elementary charge
* `full_level_structure`: OrderedDict describing properties of each energy level
  * `key::String`: Name of energy level. Spectroscopic notation encouraged, e.g. `"S1/2,f=1"`
  * `value::NamedTuple(:n, :l, :j, :f, :E)`: Quantum numbers `n`, `l`, `j`, `f`, and energy `E` (in Hz)
* `full_transitions`: Dict of all allowed transitions between energy levels
  * `key::Tuple{String,String}` Pair of levels, ordered (lower, upper) in energy
  * `value::NamedTuple(:multipole, :einsteinA)`: Leading-order multipole of the transition (e.g. `"E1"`, `"E2"`) and Einstein A coefficient (between fine structure levels only; hyperfine factors are calculated when needed)

 **Optional keywords**
 * `default_sublevel_selection`: Default value of `selected_sublevels` argument in Ion constructor
 * `gfactors`: `Dict(level::String => g::Real)` Custom Landé g-factors, if contributions from higher-than-first-order perturbations are desired
 * `nonlinear_zeeman`: `Dict` describing nonlinear contributions to Zeeman shift of certain sublevels
   * `key::Tuple{String,Real}`: sublevel name
   * `value::Function(B::Real)`: Nonlinear term(s) of Zeeman shift. Full Zeeman shift will be calculated as the sum of the usual linear term and this function
"""
Base.@kwdef struct IonProperties  # Apparently Parameters.jl is preferred for this kind of thing.
    shortname::String
    mass::Number
    charge::Number
    nuclearspin::Number
    full_level_structure::OrderedDict{
        String,
        NamedTuple{(:n, :l, :j, :f, :E), T} where T <: Tuple
    }
    full_transitions::Dict{
        Tuple{String, String},
        NamedTuple{(:multipole, :einsteinA), Tuple{String, Float64}}
    }
    default_sublevel_selection::Union{Vector{Tuple{String, String}}, Missing} = missing
    gfactors::Union{Dict{String, Number}, Missing} = missing
    nonlinear_zeeman::Union{Dict{Tuple{String, Real}, Function}, Missing} = missing
end


"""
Helper function to convert the ``full_level_structure``, as defined in the 
config, to the format expected by ``IonProperties``
"""
function process_full_level_structure(from_config::OrderedDict)::OrderedDict
    to_return = OrderedDict{String, NamedTuple{(:n, :l, :j, :f, :E)}}()
    for (level_name, level_data) in from_config
        to_return[level_name] = (
            n=level_data["n"],
            l=level_data["l"],
            j=roundnearesthalf(level_data["j"]),
            f=roundnearesthalf(level_data["f"]),
            E=level_data["E"]
        )
    end

    return to_return
end


"""
Helper function to convert the ``full_transitions``, as defined in the config,
to the format expected by ``IonProperties``.
"""
function process_full_transitions(from_config::Vector)::Dict
    full_transitions =
        Dict{Tuple{String, String}, @NamedTuple{multipole::String, einsteinA::Real}}()

    for transition_data in from_config
        key = (transition_data["from"], transition_data["to"])
        full_transitions[key] =
            (multipole=transition_data["multipole"], einsteinA=transition_data["einsteinA"])
    end

    return full_transitions
end


"""
Helper function to convert the nonlinear zeeman entries defined in the config
to the format expected by ``IonProperties``
"""
function process_nonlinear_zeeman(from_config::Vector)::Dict
    nonlinear_zeemans = Dict{Tuple{String, Integer}, Function}()
    for nonlinear_shift in from_config
        # Because of the way Julia refreshes bindings with each iteration, there
        # is no need to create closures. Wonder if this is actually efficient.
        key = (nonlinear_shift["level"], nonlinear_shift["sublevel"])
        coeffs = nonlinear_shift["coeffs"]
        nonlinear_zeemans[key] = B -> dot([B^pwr for pwr in 0:(length(coeffs)-1)], coeffs)
    end

    return nonlinear_zeemans
end
process_nonlinear_zeeman(::Missing) = missing


"""Instantiate from config file"""
function loadfromconfig(config_filepath::String)
    config = YAML.load_file(config_filepath, dicttype=OrderedDict{String, Any})

    full_level_structure = process_full_level_structure(config["full_level_structure"])
    full_transitions = process_full_transitions(config["full_transitions"])

    # Optional arguments.
    if haskey(config, "default_sublevel_selection")
        default_sublevel_selection =
            [Tuple(selection) for selection in config["default_sublevel_selection"]]
    else
        default_sublevel_selection = missing
    end
    gfactors = get(config, "gfactors", missing)
    nonlinear_zeeman = process_nonlinear_zeeman(get(config, "nonlinear_zeeman", missing))

    return IonProperties(
        shortname=config["shortname"],
        mass=config["mass"],
        charge=config["charge"],
        nuclearspin=roundnearesthalf(config["nuclearspin"]),
        full_level_structure=full_level_structure,
        full_transitions=full_transitions,
        default_sublevel_selection=default_sublevel_selection,
        gfactors=gfactors,
        nonlinear_zeeman=nonlinear_zeeman
    )
end

end  # Properties module
