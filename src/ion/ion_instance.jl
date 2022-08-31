export IonInstance

"""
    IonInstance(selected_sublevels::Vector{Tuple}[, starkshift::Dict])
Ion instance of some species

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
mutable struct IonInstance{Species <: Any} <: Ion
    # fields
    species_properties::IonProperties
    sublevels::Vector{Tuple{String, Real}}
    sublevel_aliases::Dict{String, Tuple}
    shape::Vector{Int}
    stark_shift::OrderedDict{Tuple, Real}
    ionnumber::Union{Int, Missing}
    position::Union{Real, Missing}

    # constructors (overrides default)
    function IonInstance{Species}(
        properties,
        selected_sublevels::Union{Vector{Tuple{String, T}}, String, Nothing} where {T} = nothing,
        starkshift = Dict()
    ) where {Species <: Any}
        sublevels = _construct_sublevels(selected_sublevels, properties)
        shape = [length(sublevels)]
        starkshift_full = _construct_starkshift(starkshift, sublevels)
        return new{Species}(
            properties,
            sublevels,
            Dict(),
            shape,
            starkshift_full,
            missing,
            missing
        )
    end
    function IonInstance{Species}(
        species_properties,
        sublevels,
        sublevel_aliases,
        shape,
        stark_shift,
        ionnumber,
        position
    ) where {Species <: Any}
        sublevels = deepcopy(sublevels)
        sublevel_aliases = deepcopy(sublevel_aliases)
        shape = copy(shape)
        stark_shift = deepcopy(stark_shift)
        return new{Species}(
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

function Base.print(I::IonInstance)
    println(I.species_properties.shortname + "\n")
    for (k, v) in I.selected_sublevel_structure
        println(k, ": ", v)
    end
end

Base.show(io::IO, I::IonInstance) = println(io, I.species_properties.shortname)  # suppress long output
