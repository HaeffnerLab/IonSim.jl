using .PhysicalConstants: PhysicalConstant

export Yb171

"""
    const properties_yb171
`Namedtuple` of properties of the ion species Yb171.

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
const properties_yb171 = (
    mass = 28.8384644689030595108e-26,
    charge = 1,
    nuclearspin = 1 // 2,
    full_level_structure = OrderedDict(
        # Unless otherwise indicated, frequencies/isotope shifts are as used in the Campbell group,
        #     and can be found in Conrad Roman's thesis (2020)
        # Hyperfine splitting from Fisk et al (1997)
        # Isotope shift from Martensson-Pendrill et al (1992)
        "S1/2f=0" => (n = 6, l = 0, j = 1 // 2, f = 0, E = -9.48210909315e9),
        "S1/2f=1" => (n = 6, l = 0, j = 1 // 2, f = 1, E = 3.16070303105e9),
        # Hyperfine splitting from Engelke & Tamm (1996)
        "D3/2f=0" => (n = 5, l = 2, j = 3 // 2, f = 1, E = 6.88342460964640e14),
        "D3/2f=1" => (n = 5, l = 2, j = 3 // 2, f = 2, E = 6.88350470564640e14),
        # Hyperfine splitting from Gill et al (1999)
        "D5/2f=2" => (n = 5, l = 2, j = 5 // 2, f = 2, E = 7.29477895985202e14),
        "D5/2f=3" => (n = 5, l = 2, j = 5 // 2, f = 3, E = 7.29474121985202e14),
        # Hyperfine splitting from Martensson-Pendrill et al (1992)
        # Isotope shift from Martensson-Pendrill et al (1992)
        "P1/2f=0" => (n = 6, l = 1, j = 1 // 2, f = 0, E = 8.112913749003559e14),
        "P1/2f=1" => (n = 6, l = 1, j = 1 // 2, f = 1, E = 8.112934798003559e14),
        # Hyperfine splitting from Engelke & Tamm (1996)
        "[3/2]1/2f=0" => (n = 6, l = nothing, j = 1 // 2, f = 0, E = 1.008917341058788e15),
        "[3/2]1/2f=1" => (n = 6, l = nothing, j = 1 // 2, f = 1, E = 1.008917341058788e15),
        # Hyperfine splitting can be found in Conrad Roman's thesis (2020)
        "[3/2]3/2f=1" => (n = 6, l = nothing, j = 3 // 2, f = 1, E = 8.621425511314839e14),
        "[3/2]3/2f=2" => (n = 6, l = nothing, j = 3 // 2, f = 2, E = 8.621425511314839e14),
        # Hyperfine splitting from Feldker et al (2017)
        "[5/2]5/2f=2" => (n = 6, l = nothing, j = 5 // 2, f = 2, E = 9.70461163716380e14),
        "[5/2]5/2f=3" => (n = 6, l = nothing, j = 5 // 2, f = 3, E = 9.70461163716380e14),
        # Hyperfine splitting from Gill et al (1999)
        "F7/2f=3" => (n = 6, l = 3, j = 7 // 2, f = 3, E = 6.42115934728750e14),
        "F7/2f=4" => (n = 6, l = 3, j = 7 // 2, f = 4, E = 6.42119554728750e14),

    ),
    full_transitions = Dict(
        # Einstein A Coefficients sourced from the D.R.E.A.M. database at the Univ. of Mons
        # (Quinet P., Palmeri P., Atoms 8(2), 18 (2020))
        # Laser lines
        # 411nm
        ("S1/2f=0", "D5/2f=2") => (multipole = "E2", einsteinA = 22),
        ("S1/2f=1", "D5/2f=2") => (multipole = "E2", einsteinA = 22),
        ("S1/2f=1", "D5/2f=3") => (multipole = "E2", einsteinA = 22),
        # 369nm
        ("S1/2f=0", "P1/2f=1") => (multipole = "E1", einsteinA = 1.155e8),
        ("S1/2f=1", "P1/2f=0") => (multipole = "E1", einsteinA = 1.155e8),
        ("S1/2f=1", "P1/2f=1") => (multipole = "E1", einsteinA = 1.155e8),
        # 935nm
        ("D3/2f=1", "[3/2]1/2f=0") => (multipole = "E1", einsteinA = 1.2e5),
        ("D3/2f=1", "[3/2]1/2f=1") => (multipole = "E1", einsteinA = 1.2e5),
        ("D3/2f=2", "[3/2]1/2f=1") => (multipole = "E1", einsteinA = 1.2e5),
        # 760nm
        ("F7/2f=3", "[3/2]3/2f=1") => (multipole = "E2", einsteinA = 5e4),
        ("F7/2f=3", "[3/2]3/2f=2") => (multipole = "E2", einsteinA = 5e4),
        ("F7/2f=4", "[3/2]3/2f=2") => (multipole = "E2", einsteinA = 5e4),
        # 638nm
        # 976nm
        # 861nm
        # Decay lines
        # P1/2 -> D3/2
        ("P1/2f=0", "D3/2f=1") => (multipole = "E1", einsteinA = 5.7e5),
        ("P1/2f=1", "D3/2f=2") => (multipole = "E1", einsteinA = 5.7e5),
        ("P1/2f=1", "D3/2f=2") => (multipole = "E1", einsteinA = 5.7e5),
        # [3/2]1/2 -> S1/2
        ("[3/2]1/2f=0", "S1/2f=1") => (multipole = "E1", einsteinA = 8.05e7),
        ("[3/2]1/2f=1", "S1/2f=0") => (multipole = "E1", einsteinA = 8.05e7),
        ("[3/2]1/2f=1", "S1/2f=1") => (multipole = "E1", einsteinA = 8.05e7),
        # [3/2]3/2 -> S1/2
        ("[3/2]1/2f=0", "S1/2f=1") => (multipole = "E1", einsteinA = 5.125e7),
        ("[3/2]1/2f=1", "S1/2f=0") => (multipole = "E1", einsteinA = 5.125e7),
        ("[3/2]1/2f=1", "S1/2f=1") => (multipole = "E1", einsteinA = 5.125e7),
    ),

    # Optional fields
    default_sublevel_selection = [("S1/2f=0", "all"), ("S1/2f=1", "all"), ("P1/2f=0", "all"), ("P1/2f=1", "all"), ],
    gfactors = Dict("S1/2f=0" => 1.998, "S1/2f=1" => 1.998, "D5/2f=2" => 1.202, "D5/2f=3" => 1.202, "F7/2f=3" => 1.145, "F7/2f=4" => 1.145,
    "[3/2]1/2f=0" => 1.32, "[3/2]1/2f=1" => 1.32, "[3/2]3/2f=1" => 1.44, "[3/2]3/2f=2" => 1.44,),
    nonlinear_zeeman = Dict(("S1/2f=0", 0) => B->-155.305 * B^2, ("S1/2f=1", 0) => B->155.305 * B^2),
)

"""
    Yb171(selected_sublevels::Vector{Tuple}[, starkshift::Dict])
Ion instance of species Yb171

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
mutable struct Yb171 <: Ion
    species_properties::NamedTuple
    sublevels::Vector{Tuple{String, Real}}
    sublevel_aliases::Dict{String, Tuple}
    shape::Vector{Int}
    stark_shift::OrderedDict{Tuple, Real}
    ionnumber::Union{Int, Missing}
    position::Union{Real, Missing}
    function Yb171(
        selected_sublevels::Union{Vector{Tuple{String, T}}, String, Nothing} where {T} = nothing;
        starkshift = Dict()
    )
        properties = properties_yb171

        sublevels = _construct_sublevels(selected_sublevels, properties)
        shape = [length(sublevels)]
        starkshift_full = _construct_starkshift(starkshift, sublevels)

        return new(properties, sublevels, Dict(), shape, starkshift_full, missing, missing)
    end
    # for copying
    function Yb171(
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

function Base.print(I::Yb171)
    println("171Yb\n")
    for (k, v) in I.selected_sublevel_structure
        println(k, ": ", v)
    end
end

Base.show(io::IO, I::Yb171) = println(io, "171Yb")  # suppress long output
