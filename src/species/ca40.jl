export Ca40

const properties_ca40 = (mass = PhysicalConstant(6.635943757345042e-26, "kg"),

                         charge = 1,

                         nuclearspin = 0,

                         full_level_structure = OrderedDict("S1/2" => (l=0, j=1//2, f=1//2, E=PhysicalConstant(0, "Hz")),
                                                            "D3/2" => (l=2, j=3//2, f=3//2, E=PhysicalConstant(4.09335071228e14, "Hz")),
                                                            "D5/2" => (l=2, j=5//2, f=5//2, E=PhysicalConstant(4.1115503183857306e14, "Hz")),
                                                            "P1/2" => (l=1, j=1//2, f=1//2, E=PhysicalConstant(7.554e14, "Hz")),
                                                            "P3/2" => (l=1, j=3//2, f=3//2, E=PhysicalConstant(7.621e14, "Hz")),
                                                           ),

                         default_sublevel_selection = [("S1/2", "all"),
                                                       ("D5/2", "all"),
                                                      ],

                         full_transitions = Dict(("S1/2", "D5/2") => PhysicalConstant(8.562e-1, "s⁻¹"),
                                                 ("S1/2", "P1/2") => PhysicalConstant(1.299e8, "s⁻¹"),
                                                 ("D3/2", "P1/2") => PhysicalConstant(1.060e7, "s⁻¹"),
                                                 ("S1/2", "D3/2") => PhysicalConstant(9.259e-1, "s⁻¹"),
                                                 ("S1/2", "P3/2") => PhysicalConstant(1.351e8, "s⁻¹"),
                                                 ("D3/2", "P3/2") => PhysicalConstant(1.110e6, "s⁻¹"),
                                                 ("D5/2", "P3/2") => PhysicalConstant(9.901e6, "s⁻¹"),
                                                ),

                        nonlinear_zeeman = Dict(("S1/2", -1//2) => B->1.3e-4*B^2,
                                                ("D5/2", -5//2) => B->4.5e-4*B^2,
                                               )
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
    mass::Real
    nuclearspin::Rational
    full_level_structure::OrderedDict{String,NamedTuple}
    selected_sublevel_structure::OrderedDict{Tuple,NamedTuple}
    sublevel_aliases::Dict{String,Tuple}
    shape::Vector{Int}
    full_transitions::Dict{Tuple,Real}
    selected_transitions::Dict{Tuple,Real}
    stark_shift::OrderedDict{Tuple,Real}
    nonlinear_zeeman::Dict{Tuple,Function}
    number::Union{Int,Missing}
    position::Union{Real,Missing}
    function Ca40(selected_sublevels::Union{Vector{Tuple{String,T}},String,Nothing} where T = nothing; stark_shift=Dict())
        
        properties = properties_ca40
        
        if selected_sublevels === nothing
            if :default_sublevel_selection in keys(properties)
                selected_sublevels = properties.default_sublevel_selection
            else
                @error "no level structure specified in constructor, and no default level structure specified for this ion species"
            end
        elseif selected_sublevels == "all"
            selected_sublevels = [(sublevel, "all") for sublevel in keys(properties.full_level_structure)]
        end
        m = properties.mass
        i = properties.nuclearspin
        fls = properties.full_level_structure
        ft = properties.full_transitions
        sss, st = _structure(selected_sublevels, fls, ft)
        shape = [length(sss)]
        ss_full = OrderedDict{Tuple,Real}()
        nlz = OrderedDict{Tuple,Function}()
        for sublevel in keys(sss)
            haskey(stark_shift, sublevel) ? ss_full[sublevel] = stark_shift[sublevel] : ss_full[sublevel] = 0.
            haskey(properties.nonlinear_zeeman, sublevel) ? nlz[sublevel] = properties.nonlinear_zeeman[sublevel] : nlz[sublevel] = B->0.
        end
        new(m, i, fls, sss, Dict(), shape, ft, st, ss_full, nlz, missing, missing)
    end
    # for copying
    function Ca40(  
            mass, nuclearspin, full_level_structure, selected_sublevel_structure, sublevel_aliases, shape,
            full_transitions, selected_transitions, stark_shift, nonlinear_zeeman, number, position
        )
        selected_sublevel_structure = deepcopy(selected_sublevel_structure)
        sublevel_aliases = deepcopy(sublevel_aliases)
        shape = copy(shape)
        selected_transitions = deepcopy(selected_transitions)
        stark_shift = deepcopy(stark_shift)
        nonlinear_zeeman = deepcopy(nonlinear_zeeman)
        new(mass, nuclearspin, full_level_structure, selected_sublevel_structure, sublevel_aliases, shape,
        full_transitions, selected_transitions, stark_shift, nonlinear_zeeman, number, position)
    end
end