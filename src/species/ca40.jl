export Ca40

const properties_ca40 = (mass = PhysicalConstant(6.635943757345042e-26, "kg"),

                         charge = 1,

                         i = 0,

                         full_level_structure = OrderedDict(
                                                            "S1/2" => (l=0, j=1//2, f=1//2, E=PhysicalConstant(0, "Hz")),
                                                            "D3/2" => (l=2, j=3//2, f=3//2, E=PhysicalConstant(4.09335071228e14, "Hz")),
                                                            "D5/2" => (l=2, j=5//2, f=5//2, E=PhysicalConstant(4.1115503183857306e14, "Hz")),
                                                            "P1/2" => (l=1, j=1//2, f=1//2, E=PhysicalConstant(7.554e14, "Hz")),
                                                            "P3/2" => (l=1, j=3//2, f=3//2, E=PhysicalConstant(7.621e14, "Hz")),
                                                            ),

                         full_transitions = Dict(
                                                 ("S1/2", "D5/2") => PhysicalConstant(8.562e-1, "s⁻¹"),
                                                 ("S1/2", "P1/2") => PhysicalConstant(1.299e8, "s⁻¹"),
                                                 ("D3/2", "P1/2") => PhysicalConstant(1.060e7, "s⁻¹"),
                                                 ("S1/2", "D3/2") => PhysicalConstant(9.259e-1, "s⁻¹"),
                                                 ("S1/2", "P3/2") => PhysicalConstant(1.351e8, "s⁻¹"),
                                                 ("D3/2", "P3/2") => PhysicalConstant(1.110e6, "s⁻¹"),
                                                 ("D5/2", "P3/2") => PhysicalConstant(9.901e6, "s⁻¹"),
                                                ),

                        nonlin_zeeman = Dict(
                                             ("S1/2", -1//2) => B->1.3e-4*B^2,
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
    level_structure::OrderedDict{String,NamedTuple}
    selected_level_structure::OrderedDict{String,NamedTuple}
    shape::Vector{Int}
    matrix_elements::OrderedDict{Tuple,Function}
    selected_matrix_elements::OrderedDict{Tuple,Function}
    stark_shift::OrderedDict{String,Real}
    number::Union{Int,Missing}
    position::Union{Real,Missing}
    function Ca40(selected_level_structure; ss=Dict())
        fls, sls_dict, me, me_dict=_structure(selected_level_structure)
        shape = [length(sls_dict)]
        ss_full = OrderedDict{String,Float64}()
        for level in keys(sls_dict)
            haskey(ss, level) ? ss_full[level] = ss[level] : ss_full[level] = 0.
        end
        new(m_ca40, fls, sls_dict, shape, me, me_dict, ss_full, missing, missing)
    end
    Ca40(;ss=Dict()) = Ca40("default", ss=ss)
    # for copying
    function Ca40(  
            mass, level_structure, selected_level_structure, shape, matrix_elements,
            selected_matrix_elements, stark_shift, number, position
        )
        selected_level_structure = deepcopy(selected_level_structure)
        shape = copy(shape)
        selected_matrix_elements = deepcopy(selected_matrix_elements)
        stark_shift = deepcopy(stark_shift)
        new(mass, level_structure, selected_level_structure, shape, matrix_elements, 
            selected_matrix_elements, stark_shift, number, position)
    end
end