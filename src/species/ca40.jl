using .PhysicalConstants:PhysicalConstant

export Ca40

const properties_ca40 = (mass = PhysicalConstant(6.635943757345042e-26, "kg"),

                         charge = 1,

                         nuclearspin = 0,

                         full_level_structure = OrderedDict("S1/2" => (n=4, l=0, j=1//2, f=1//2, E=0),
                                                            "D3/2" => (n=3, l=2, j=3//2, f=3//2, E=4.09335071228e14),
                                                            "D5/2" => (n=3, l=2, j=5//2, f=5//2, E=4.1115503183857306e14),
                                                            "P1/2" => (n=4, l=1, j=1//2, f=1//2, E=7.554e14),
                                                            "P3/2" => (n=4, l=1, j=3//2, f=3//2, E=7.621e14),
                                                           ),

                         full_transitions = Dict(("S1/2", "D5/2") => (multipole="E2", einsteinA=8.562e-1),
                                                 ("S1/2", "P1/2") => (multipole="E1", einsteinA=1.299e8),
                                                 ("D3/2", "P1/2") => (multipole="E1", einsteinA=1.060e7),
                                                 ("S1/2", "D3/2") => (multipole="E2", einsteinA=9.259e-1),
                                                 ("S1/2", "P3/2") => (multipole="E1", einsteinA=1.351e8),
                                                 ("D3/2", "P3/2") => (multipole="E1", einsteinA=1.110e6),
                                                 ("D5/2", "P3/2") => (multipole="E1", einsteinA=9.901e6),
                                                ),
                        
                         # Optional fields
                         default_sublevel_selection = [("S1/2", "all"),
                                                       ("D5/2", "all"),
                                                      ],

                         gfactors = Dict("S1/2" => 2.00225664,
                                         "D5/2" => 1.2003340
                                        ),

                         #nonlinear_zeeman = Dict(("S1/2", -1//2) => B->1.3e-4*B^2,
                         #                        ("D5/2", -5//2) => B->4.5e-4*B^2,
                         #                       )
                         )







"""
    Ca40(selected_level_structure::Vector{String}[, stark_shift])

#### user-defined fields
* `selected_level_structure`:
    keys ⊂ `["S-1/2", "S+1/2", "D-5/2", "D-3/2", "D-1/2", "D+1/2", "D+3/2", "D+5/2"]`.
    Values are a `NamedTuple` with:
    * `l`: orbital angular momentum
    * `j`: total angular momentum
    * `mⱼ`: projection of total angular momentum along quantization axis
    * `E`: relative energies
    Note: indexing the instantiated structure with one of these strings will return
    the corresponding `Ket`.
* `stark_shift`: A dictionary with keys denoting the selected levels and values, a real
    number for describing a shift of the level's energy. This is just a convenient way to add
    Stark shifts to the simulation without additional resources.
#### fixed fields
* `mass::Real`: The ion's mass in kg.
* `level_structure`: A full description of the ion's electronic structure.
* `matrix_elements::OrderedDict{Tuple,Function}`: Same as `selected_matrix_elements` but for
    all of the ion's allowed transitions.
#### derived fields
* `selected_matrix_elements`: Functions for the allowed transitions (contained in the
    selected levels) that return the corresponding coupling strengths. These functions take
    as arguments:
    * `Efield`: magnitude of the electric field at the position of the ion [V/m]
    * `γ`: ``ϵ̂⋅B̂`` (angle between laser polarization and B-field)
    * `ϕ`: ``k̂⋅B̂`` (angle between laser k-vector and B-field)
* `shape::Vector{Int}`: Indicates the dimension of the used Hilbert space.
* `number`: When the ion is added to an `IonConfiguration`, this value keeps track of its
    order in the chain.
* `position`: @hen the ion is added to an `IonConfiguration`, this value keeps track of its
    physical position in meters.
"""
mutable struct Ca40 <: Ion
    species_properties::NamedTuple
    sublevels::Vector{Tuple{String,Real}}
    sublevel_aliases::Dict{String,Tuple}
    shape::Vector{Int}
    stark_shift::OrderedDict{Tuple,Real}
    ionnumber::Union{Int,Missing}
    position::Union{Real,Missing}
    function Ca40(selected_sublevels::Union{Vector{Tuple{String,T}},String,Nothing} where T=nothing; starkshift=Dict())
        
        properties = properties_ca40
        
        sublevels = _construct_sublevels(selected_sublevels, properties)
        shape = [length(sublevels)]
        starkshift_full = _construct_starkshift(starkshift, sublevels)

        new(properties, sublevels, Dict(), shape, starkshift_full, missing, missing)
    end
    # for copying
    function Ca40(species_properties, sublevels, sublevel_aliases, shape, stark_shift, ionnumber, position)
        sublevels = deepcopy(sublevels)
        sublevel_aliases = deepcopy(sublevel_aliases)
        shape = copy(shape)
        stark_shift = deepcopy(stark_shift)
        new(species_properties, sublevels, sublevel_aliases, shape, stark_shift, ionnumber, position)
    end
end