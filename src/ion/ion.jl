using IonSim.PhysicalConstants: e

export Ion

"""
    Ion
The physical parameters defining an isolated ion's internal structure.
"""
abstract type Ion <: IonSimBasis end

# Expected fields
speciesproperties(I::Ion)::IonProperties = I.species_properties
sublevels(I::Ion)::Vector{Tuple{String, Real}} = I.sublevels
sublevel_aliases(I::Ion)::Dict{String, Tuple} = I.sublevel_aliases
shape(I::Ion)::Vector{Int} = I.shape
stark_shift(I::Ion)::OrderedDict{Tuple, Real} = I.stark_shift
ionnumber(I::Ion)::Union{Int, Missing} = I.ionnumber
ionposition(I::Ion)::Union{Real, Missing} = I.position

# Properties of Ion species
mass(I::Ion)::Real = speciesproperties(I).mass
charge(I::Ion)::Real = speciesproperties(I).charge * e
nuclearspin(I::Ion)::Rational = speciesproperties(I).nuclearspin

# Handle Ion state differently
Base.getindex(I::Ion, state::Union{Tuple{String, Real}, String, Int}) = ionstate(I, state)

# Warn if Ion not in configuration
function Base.getproperty(I::Ion, s::Symbol)
    if s == :ionnumber || s == :position
        if typeof(getfield(I, s)) <: Missing
            @warn "ion has not been added to a configuration"
            return missing
        end
    end
    return getfield(I, s)
end

# Ions are equal when they are the same species, contain the same sublevels, and have the same stark shifts
# Sublevel aliases don't matter; nor does sublevel ordering
function Base.:(==)(b1::T, b2::T) where {T <: Ion}
    return (
        b1.species_properties == b2.species_properties &&
        sort(b1.sublevels) == sort(b2.sublevels) &&
        sort(b1.stark_shift) == sort(b2.stark_shift)
    )
end

# stark shift constructor
function _construct_starkshift(starkshift, sublevels)
    starkshift_full = OrderedDict{Tuple, Real}()
    for sublevel in sublevels
        starkshift_full[sublevel] =
            (haskey(starkshift, sublevel) ? starkshift[sublevel] : 0.0)
    end
    return starkshift_full
end

# sublevels constructor
function _construct_sublevels(selected_sublevels, properties)
    full_level_structure = properties.full_level_structure

    # If selected_sublevels is blank, use the default selection. If it is "all", use all sublevels.
    if selected_sublevels === nothing
        if !ismissing(properties.default_sublevel_selection)
            selected_sublevels = properties.default_sublevel_selection
        else
            @error "no level structure specified in constructor, and no default level structure specified for this ion species"
        end
    elseif selected_sublevels == "all"
        selected_sublevels = [(sublevel, "all") for sublevel in keys(full_level_structure)]
    end

    # Construct the list of sublevels
    sublevels = []
    for manifold in selected_sublevels
        # Ensure that the string is a valid level
        level = manifold[1]
        @assert level in keys(full_level_structure) "invalid level $level"
        @assert level ∉ [k[1] for k in keys(sublevels)] "multiple instances of level $level in ion constructor call"
        level_structure = full_level_structure[level]

        # Add chosen sublevels
        selectedms = manifold[2]
        f = level_structure.f

        m_allowed = Array((-f):f)
        if selectedms == "all"
            selectedms = m_allowed
        elseif !(typeof(selectedms) <: Array)
            selectedms = [selectedms]
        end
        for m in selectedms
            m = Rational(m)
            @assert m in m_allowed "Zeeman sublevel m = $m not allowed for state $level with f = $f"
            @assert (level, m) ∉ sublevels "repeated instance of sublevel $m in state $level"
            push!(sublevels, (level, m))
        end
    end

    return sublevels
end

# Update Ion sublevels properly
function Base.setproperty!(I::Ion, s::Symbol, v::Tv) where {Tv}
    if (s == :species_properties || s == :shape || s == :number || s == :position)
        return
    elseif s == :sublevels
        Core.setproperty!(I, :sublevels, _construct_sublevels(v, speciesproperties(I)))
        # Also update the stark shift dict as necessary; keep old stark shift values and assign zero stark shift to new sublevels
        starkshift_full_old = stark_shift(I)
        starkshift_full_new = OrderedDict{Tuple, Real}()
        for sublevel in sublevels(I)
            starkshift_full_new[sublevel] = (
                haskey(starkshift_full_old, sublevel) ? starkshift_full_old[sublevel] : 0.0
            )
        end
        Core.setproperty!(I, :stark_shift, starkshift_full_new)
        return
    elseif s == :sublevel_aliases
        for (alias, sublevel) in v
            validatesublevel(I, sublevel)
        end
        Core.setproperty!(I, :sublevel_aliases, v)
    elseif s == :stark_shift
        starkshift_full_new = stark_shift(I) # New stark shift dict is initially identical to the old one
        for (sublevel, shift) in v
            validatesublevel(I, sublevel)
            starkshift_full_new[sublevel] = shift # Change the shifts as necessary
        end
        Core.setproperty!(I, :stark_shift, starkshift_full_new)
    end
end
