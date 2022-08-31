export IonProperties

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
struct IonProperties
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
    default_sublevel_selection::Union{Vector{Tuple{String, String}}, Missing}
    gfactors::Union{Dict{String, Number}, Missing}
    nonlinear_zeeman::Union{Dict{Tuple{String, Real}, Function}, Missing}
end
IonProperties(;
    shortname,
    mass,
    charge,
    nuclearspin,
    full_level_structure,
    full_transitions,
    default_sublevel_selection = missing,
    gfactors = missing,
    nonlinear_zeeman = missing
) = IonProperties(
    shortname,
    mass,
    charge,
    nuclearspin,
    full_level_structure,
    full_transitions,
    default_sublevel_selection,
    gfactors,
    nonlinear_zeeman
)
