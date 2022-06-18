export Be9

"""
    const properties_be9
`Namedtuple` of properties of the ion species Be9.

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
const properties_be9 = (
    mass = 1.496508080073e-26,
    charge = 1,
    nuclearspin = 3 // 2,
    full_level_structure = OrderedDict(
        "S1/2f=1" => (n = 2, l = 0, j = 1 // 2, f = 1, E = 0.78126104631e9),
        "S1/2f=2" => (n = 2, l = 0, j = 1 // 2, f = 2, E = -0.468756627786e9),
        "P1/2f=1" => (n = 2, l = 1, j = 1 // 2, f = 1, E = 957.4772214787497e12),
        "P1/2f=2" => (n = 2, l = 1, j = 1 // 2, f = 2, E = 957.4769842787498e12),
        "P3/2f=0" => (n = 2, l = 1, j = 3 // 2, f = 0, E = 957.6770501980199e12),
        "P3/2f=1" => (n = 2, l = 1, j = 3 // 2, f = 1, E = 957.6770491780198e12),
        "P3/2f=2" => (n = 2, l = 1, j = 3 // 2, f = 2, E = 957.6770471380199e12),
        "P3/2f=3" => (n = 2, l = 1, j = 3 // 2, f = 3, E = 957.6770440780199e12),
    ),
    full_transitions = Dict(
        ("S1/2f=1", "P1/2f=1") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=1", "P1/2f=2") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=2", "P1/2f=1") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=2", "P1/2f=2") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=1", "P3/2f=0") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=1", "P3/2f=1") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=1", "P3/2f=2") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=2", "P1/2f=1") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=2", "P1/2f=2") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=2", "P1/2f=3") => (multipole = "E1", einsteinA = 19.4e6),
    ),

    # Optional fields
    default_sublevel_selection = [
        ("S1/2f=1", "all"),
        ("S1/2f=2", "all"),
        ("P1/2f=1", "all"),
        ("P1/2f=2", "all"),
        ("P3/2f=0", "all"),
        ("P3/2f=1", "all"),
        ("P3/2f=2", "all"),
        ("P3/2f=3", "all"),
    ],
)

"""
    Be9(selected_sublevels::Vector{Tuple}[, starkshift::Dict])
Ion instance of species Be9

`selected_sublevels` specifies which energy sublevels will be present in the Hilbert space of this Ion instance, as a subset of all possible sublevels.

Each element of `selected_sublevels` is a 2-element Tuple (level, sublevels), with the first element being the name of a level and the second specifying which sublevels should be included.
Allowed sublevels are those whose magnetic quantum number `m` is in the set {`-f`, `-f+1`, `-f+2`, ... `f-1`, `f`}, where `f` is the total angular momentum quantum number of `level`.
For each `level` specified there are three allowed options to specify the set of `sublevels` to include:
* `sublevels::Real`: Includes only one `m = sublevels`
* `sublevels::Vector{Real}`: Includes all sublevels whose magnetic quantum number `m` is in `sublevels`
* `sublevels = "all"`: Includes all allowed sublevels

If instead `selected_sublevels = "all"`, then all sublevels of all levels are included.

Omission of a level in `selected_sublevels` will exclude all sublevels.

**Fields**
* `species_properties::NamedTuple`: Contains constants specifying parameters specific to species
* `sublevels`::Vector{Tuple{String,Real}}: List of all sublevels present in the Hilbert space
* `sublevel_aliases::Dict{String,Tuple}`: Dict specifying aliases assigned to sublevels, in the format `alias => sublevel`
* `shape`::Vector{Int}: Dimension of the Hilbert space
* `stark_shift::OrderedDict`: A dictionary with keys denoting the selected levels and values, a real number for describing a shift of the level's energy. This is just a convenient way to add Stark shifts to the simulation without additional resources
* `ionnumber`: When the ion is added to an `IonConfiguration`, this value keeps track of its order in the chain
* `position`: When the ion is added to an `IonConfiguration`, this value keeps track of its physical position in meters
"""
mutable struct Be9 <: Ion
    species_properties::NamedTuple
    sublevels::Vector{Tuple{String, Real}}
    sublevel_aliases::Dict{String, Tuple}
    shape::Vector{Int}
    stark_shift::OrderedDict{Tuple, Real}
    ionnumber::Union{Int, Missing}
    position::Union{Real, Missing}
    function Be9(
        selected_sublevels::Union{Vector{Tuple{String, T}}, String, Nothing} where {T} = nothing;
        starkshift = Dict()
    )
        properties = properties_be9

        sublevels = _construct_sublevels(selected_sublevels, properties)
        shape = [length(sublevels)]
        starkshift_full = _construct_starkshift(starkshift, sublevels)

        return new(properties, sublevels, Dict(), shape, starkshift_full, missing, missing)
    end
    # for copying
    function Be9(
        species_properties,
        sublevels,
        sublevel_aliases,
        shape,
        stark_shift,
        ionnumber,
        position
    )
        sublevels = deepcopy(sublevels)
        sublevel_aliases = deepcopy(sublevel_aliases)
        shape = copy(shape)
        stark_shift = deepcopy(stark_shift)
        return new(
            species_properties,
            sublevels,
            sublevel_aliases,
            shape,
            stark_shift,
            ionnumber,
            position
        )
    end
end

function Base.print(I::Be9)
    println("⁹Be\n")
    for (k, v) in I.selected_sublevel_structure
        println(k, ": ", v)
    end
end

Base.show(io::IO, I::Be9) = println(io, "⁹Be")  # suppress long output
