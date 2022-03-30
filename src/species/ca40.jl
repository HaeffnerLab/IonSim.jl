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
This needs a docstring
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