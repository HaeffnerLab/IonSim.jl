using Unitful
using .PhysicalConstants: INVERSE_TIME
export Ca40

"""
    const properties_ca40
`Namedtuple` of properties of the ion species Ca40.

**Required keywords**
* `mass`: Mass of ion
* `charge`: Charge of ion
* `full_level_structure`: OrderedDict describing properties of each energy level
  * `key::String`: Name of energy level. Spectroscopic notation encouraged, e.g. `"S1/2,f=1"`
  * `value::NamedTuple(:n, :l, :j, :f, :E)`: Quantum numbers `n`, `l`, `j`, `f`, and energy `E`
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
const properties_ca40 = (
    mass = 6.635943757345042e-26u"kg",
    charge = 1u"q",
    nuclearspin = 0,
    full_level_structure = OrderedDict(
        "S1/2" => (n = 4, l = 0, j = (1 // 2), f = (1 // 2), E = 0u"1/s"),
        "D3/2" => (n = 3, l = 2, j = (3 // 2), f = (3 // 2), E = 4.09335071228e14u"1/s"),
        "D5/2" =>
            (n = 3, l = 2, j = (5 // 2), f = (5 // 2), E = 4.1115503183857306e14u"1/s"),
        "P1/2" => (n = 4, l = 1, j = (1 // 2), f = (1 // 2), E = 7.554e14u"1/s"),
        "P3/2" => (n = 4, l = 1, j = (3 // 2), f = (3 // 2), E = 7.621e14u"1/s"),
    ),
    full_transitions = Dict(
        ("S1/2", "D5/2") => (multipole = "E2", einsteinA = 8.562e-1u"1/s"),
        ("S1/2", "P1/2") => (multipole = "E1", einsteinA = 1.299e8u"1/s"),
        ("D3/2", "P1/2") => (multipole = "E1", einsteinA = 1.060e7u"1/s"),
        ("S1/2", "D3/2") => (multipole = "E2", einsteinA = 9.259e-1u"1/s"),
        ("S1/2", "P3/2") => (multipole = "E1", einsteinA = 1.351e8u"1/s"),
        ("D3/2", "P3/2") => (multipole = "E1", einsteinA = 1.110e6u"1/s"),
        ("D5/2", "P3/2") => (multipole = "E1", einsteinA = 9.901e6u"1/s"),
    ),

    # Optional fields
    default_sublevel_selection = [("S1/2", "all"), ("D5/2", "all"),],
    gfactors = Dict("S1/2" => 2.00225664, "D5/2" => 1.2003340),

    #nonlinear_zeeman = Dict(("S1/2", -1//2) => B->1.3e-4*B^2,
    #                        ("D5/2", -5//2) => B->4.5e-4*B^2)  # Syntax example, not numerically accurate
)

"""
    Ca40(selected_sublevels::Vector{Tuple}[, starkshift::Dict])
Ion instance of species Ca40

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
* `position`: When the ion is added to an `IonConfiguration`, this value keeps track of its physical position
"""
mutable struct Ca40 <: Ion
    species_properties::NamedTuple
    sublevels::Vector{Tuple{String, Real}}
    sublevel_aliases::Dict{String, Tuple}
    shape::Vector{Int}
    stark_shift::OrderedDict{Tuple, INVERSE_TIME}
    ionnumber::Union{Int, Missing}
    position::Union{Unitful.Length, Missing}
    function Ca40(
        selected_sublevels::Union{Vector{Tuple{String, T}}, String, Nothing} where {T} = nothing;
        starkshift = Dict()
    )
        properties = properties_ca40

        sublevels = _construct_sublevels(selected_sublevels, properties)
        shape = [length(sublevels)]
        starkshift_full = _construct_starkshift(starkshift, sublevels)

        return new(properties, sublevels, Dict(), shape, starkshift_full, missing, missing)
    end
    # for copying
    function Ca40(
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

function Base.print(I::Ca40)
    println("⁴⁰Ca\n")
    for (k, v) in I.selected_sublevel_structure
        println(k, ": ", v)
    end
end

Base.show(io::IO, I::Ca40) = println(io, "⁴⁰Ca")  # suppress long output
