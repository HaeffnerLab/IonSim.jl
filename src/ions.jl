using WignerSymbols: wigner3j
using .PhysicalConstants: e, ca40_qubit_transition_frequency, m_ca40, ħ, α, μB


export mass, nuclearspin, full_level_structure, selected_sublevel_structure, sublevel_aliases,
       full_transitions, selected_transitions, get_basis, stark_shift, ionnumber, ion_position,
       set_sublevel_alias!, gJ, zeeman_shift, matrix_elements, zero_stark_shift, Ion




#############################################################################################
# Ion - the physical parameters defining an ion's structure
#############################################################################################

"""
    Ion
The physical parameters defining an isolated ion's internal structure.
"""
abstract type Ion <: IonSimBasis end




#############################################################################################
# Object fields
#############################################################################################

speciesproperties(I::Ion)::NamedTuple = I.species_properties
ionsublevels(I::Ion)::OrderedDict{Tuple,NamedTuple} = I.sublevels
sublevel_aliases(I::Ion)::Dict{String,Tuple} = I.sublevel_aliases
shape(I::Ion)::Vector{Int} = I.shape
stark_shift(I::Ion)::OrderedDict{Tuple,Real} = I.stark_shift
ionnumber(I::Ion)::Union{Int,Missing} = I.number
ionposition(I::Ion)::Union{Real,Missing} = I.position




#############################################################################################
# General properties of species
#############################################################################################

mass(I::Ion)::Real = speciesproperties(I).mass
charge(I::Ion)::Real = speciesproperties(I).charge * e
nuclearspin(I::Ion)::Rational = speciesproperties(I).nuclearspin




#############################################################################################
# Functions to modify ion properties
#############################################################################################

"""
This needs a docstring
"""
function zero_stark_shift!(I::Ion)
    for sublevel in keys(stark_shift(I))
        I.stark_shift[sublevel] = 0.0
    end
end


validatesublevel(I::Ion, sublevel::Tuple{String,Real}) = @assert sublevel in sublevels(I) "ion does not contain sublevel $sublevel"
validatesublevel(I::Ion, alias::String) = validatesublevel(I, alias2sublevel(I, alias))


"""
This needs a docstring
"""
function set_sublevel_alias!(I::Ion, sublevel::Tuple{String,Real}, alias::String)
    validatesublevel(I, sublevel)
    @assert alias ∉ levels(I) "cannot make alias name identical to level name ($alias)"
    I.sublevel_aliases[alias] = sublevel
end
function set_sublevel_alias!(I::Ion, pairs::Vector{Tuple{Tuple{String,Real},String}})
    for (sublevel, alias) in pairs
        set_sublevel_alias!(I, sublevel, alias)
    end
end


function alias2sublevel(I::Ion, alias::String)
    all_aliases = I.sublevel_aliases
    @assert alias in keys(all_aliases) "no sublevel with alias $alias"
    return all_aliases[alias]
end



#############################################################################################
# Properties of ion electronic levels and sublevels
#############################################################################################

"""
This needs a docstring
"""
ionlevels(I::Ion) = unique([sublevel[1] for sublevel in sublevels(I)])


"""
This needs a docstring
"""
function quantumnumbers(I::Ion, sublevel::Tuple{String,Real})
    validatesublevel(I, sublevel)
    levelstruct = speciesproperties(I).full_level_structure[level]
    names = (:n, :i, :s, :l, :j, :f, :m)
    values = [levelstruct.n, nuclearspin(I), 1//2, levelstruct.l, levelstruct.j, levelstruct.f, Rational(sublevel[2])]
    return NamedTuple{names}[values]
end
function quantumnumbers(I::Ion, level_or_alias::String)
    # If the second argument is a String, it could be either a level name or the alias of a sublevel
    if level_or_alias in levels(I)
        # Second argument is a level name. Leave out the m quantum number
        levelstruct = speciesproperties(I).full_level_structure[level]
        names = (:n, :i, :s, :l, :j, :f)
        values = [levelstruct.n, nuclearspin(I), 1//2, levelstruct.l, levelstruct.j, levelstruct.f]
        return NamedTuple{names}[values]
    else
        # Second argument is a sublevel alias.
        quantumnumbers(I, alias2sublevel(I, level_or_alias))
    end
end


"""
    landegf(I::Ion, level::String)
Landé g-factor

### args
###############################FILL THIS IN######################
"""
landegf(l::Real, j::Real, f::Real, i::Real, s::Real=1//2) = landegj(l, j, s)/2 * (1 + ((j*(j+1) - i*(i+1)) / (f*(f+1))))
landegf(qnums::NamedTuple) = landegf(qnums.l, qnums.j, qnums.f, qnums.i, qnums.s)
function landegf(I::Ion, level::String)
    properties = speciesproperties(I)
    if haskey(properties, :gfactors) && haskey(properties.gfactors, level)
        return properties.gfactors[level]
    else
        return landegf(quantumnumbers(I, level))
    end
end
landegj(l::Real, j::Real, s::Real=1//2) = 3//2 + (s*(s+1) - l*(l+1)) / (2j*(j+1))


"""
This needs a docstring
(Don't forget that there's a method above for just stark_shift(I::Ion))
"""
stark_shift(I::Ion, sublevel::Tuple{String,Real}) = stark_shift(I)[sublevel]
stark_shift(I::Ion, alias::String) = stark_shift(I, alias2sublevel(I, alias))


"""
NEEDS TO BE CHANGED
    zeeman_shift(B::Real, l::Real, j::Real, mⱼ::Real)
``ΔE = (μ_B/ħ) ⋅ g_J(l, j) ⋅ B ⋅ mⱼ / 2π``
### args
* `B`: magnitude of B-field at ion
* `l`: orbital angular momentum quantum number
* `j`: total angular momentum quantum number
* `mⱼ`: projection of total angular momentum along quantization axis
"""
zeeman_shift(B::Real, g::Real, m::Real) = (μB/ħ) * g * B * m / 2π
zeeman_shift(B::Real, l::Real, j::Real, f::Real, m::Real, i::Real, s::Real=1//2) = zeeman_shift(B, landegf(l, j, f, i, s), m)
zeeman_shift(B::Real, qnums::NamedTuple) = zeeman_shift(B, qnums.l, qnums.j, qnums.f, qnums.m, qnums.i, qnums.s)
function zeeman_shift(I::Ion, sublevel::Union{Tuple{String,Real},String}, B::Real)
    validatesublevel(I, sublevel)
    properties = speciesproperties(I)
    if haskey(properties, :nonlinear_zeeman) && haskey(properties.nonlinear_zeeman, sublevel)
        nonlinear = properties.nonlinear_zeeman[sublevel](B)
    else
        nonlinear = 0.0
    return zeeman_shift(B, landegf(I, sublevel[1]), sublevel[2]) + nonlinear
end
zeeman_shift(I::Ion, alias::String, B::Real) = zeeman_shift(I, alias2sublevel(I, alias), B)


"""
This needs a docstring
"""
function energy(I::Ion, sublevel::Tuple{String,Real}; B=0, ignore_starkshift=false)
    validatesublevel(I, sublevel)
    E0 = speciesproperties(I).full_level_structure[sublevel[0]].E
    zeeman = zeeman_shift(I, sublevel, B)
    stark = (ignore_starkshift ? 0.0 : stark_shift(I, sublevel))
    return E0 + zeeman + stark
end
function energy(I::Ion, level_or_alias::String; B=0, ignore_starkshift=false)
    # If the second argument is a String, it could be either a level name or the alias of a sublevel
    if level_or_alias in levels(I)
        # Second argument is a level name. Return the bare energy of that level.
        return speciesproperties(I).full_level_structure[level_or_alias].E
    else
        # Second argument is a sublevel alias.
        return energy(I, alias2sublevel(I, level_or_alias), B=B, ignore_starkshift=ignore_starkshift)
    end
end


"""
This needs a docstring
"""
function transitionfrequency(I::Ion, sublevel1::Tuple{String,Real}, sublevel2::Tuple{String,Real}, B::Real; ignore_starkshift=false)
    return abs(energy(I, sublevel1, B=B, ignore_starkshift=ignore_starkshift) - energy(I, sublevel2, B=B, ignore_starkshift=ignore_starkshift))
end


"""
This needs a docstring
"""
function ionleveltransitions(I::Ion)
    list = []
    levels = ionlevels(I)
    for levelpair in keys(speciesproperties(I).full_transitions)
        if levelpair[1] in levels && levelpair[2] in levels
            push!(list, levelpair)
        end
    end
    return list
end


"""
This needs a docstring
"""
function ionsubleveltransitions(I::Ion)
    list = []
    leveltransitions = ionleveltransitions(I)
    for transition in leveltransitions
        (L1, L2) = transition
        sublevels1 = [sublevel for sublevel in ionsublevels(I) if sublevel[0]==L1]
        sublevels2 = [sublevel for sublevel in ionsublevels(I) if sublevel[0]==L2]
        for sl1 in sublevels1
            for sl2 in sublevels2
                push!(list, (sl1, sl2))
            end
        end
    end
    return list
end


"""
This needs a docstring
"""
function einsteinA(I::Ion, L1::Tuple{String,Real}, L2::Tuple{String,Real})
    @assert (L1, L2) in ionleveltransitions(I) "invalid transition $L1 -> $L2"
    return speciesproperties(I).full_transitions[(L1, L2)]
end
einsteinA(I::Ion, L1::Tuple{String,Real}, L2::String) = einsteinA(I, L1, alias2sublevel(I, L2))
einsteinA(I::Ion, L1::String, L2::Tuple{String,Real}) = einsteinA(I, (alias2sublevel(I, L1), L2))
einsteinA(I::Ion, L1::String, L2::String) = einsteinA(I, (alias2sublevel(I, L1), alias2sublevel(I, L2)))
einsteinA(I::Ion, Lpair::Tuple) = einsteinA(I, Lpair[1], Lpair[2])


"""
This needs a docstring
Main method is currently a placeholder
"""
function matrix_element(Δl::Int, j1::Real, j2::Real, f1::Real, f2::Real, Δm::Int, ΔE::Real, A12::Real, Efield::Function, khat::NamedTuple, ϵhat::NamedTuple, Bhat::NamedTuple=(;z=1))
    # Decide type of transition
    # Rotate unit vectors so that B is in z-direction?
    # Calculate matrix element
    # Return a function of time
end
function matrix_element(I::Ion, transition::Tuple, Efield::Function, khat::NamedTuple, ϵhat::NamedTuple, Bhat::NamedTuple=(;z=1))
    qn1 = quantumnumbers(I, transition[1])
    qn2 = quantumnumbers(I, transition[2])
    E1 = energy(I, transition[1], ignore_starkshift=true)
    E2 = energy(I, transition[2], ignore_starkshift=true)
    matrix_element(qn2.l-qn1.l, qn1.j, qn2.j, qn1.f, qn2.f, qn2.m-qn1.m, abs(E2-E1), einsteinA(I, transition), Efield, khat, ϵhat, Bhat)
end
matrix_element(I::Ion, transition::Tuple, T::Trap, laser::Laser) = matrix_element(I, transition, laser.E, laser.k, laser.ϵ, T.Bhat)




#############################################################################################
# Functions for constructing ion objects
#############################################################################################

function _construct_sublevels(selected_sublevels, properties)
    full_level_structure = properties.full_level_structure

    # If selected_sublevels is blank, use the default selection. If it is "all", use all sublevels.
    if selected_sublevels === nothing
        if :default_sublevel_selection in keys(properties)
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
        @assert level ∉ [k[1] for k in keys(sublevel_structure)] "multiple instances of level $level in ion constructor call"
        level_structure = full_level_structure[level]

        # Add chosen sublevels
        selectedms = manifold[2]
        f = level_structure.f

        m_allowed = Array(-f:f)
        if selectedms == "all"
            selectedms = m_allowed
        elseif !(typeof(selectedms) <: Array)
            selectedms = [selectedms]
        end
        for m in selectedms
            m = Rational(m)
            @assert m in m_allowed "Zeeman sublevel m = $m not allowed for state $level with f = $f"
            @assert (level, m) ∉ keys(sublevel_structure) "repeated instance of sublevel $m in state $level"
            push!(sublevels, (level, m))
        end
    end

    return sublevels
end

function _construct_starkshift(starkshift, sublevels)
    starkshift_full = OrderedDict{Tuple,Real}()
    for sublevel in sublevels
        starkshift_full[sublevel] = (haskey(starkshift, sublevel) ? starkshift[sublevel] : 0.)
    end
    return starkshift_full
end



#############################################################################################
# Overrides of Base functions
#############################################################################################

function Base.print(I::Ca40)
    println("⁴⁰Ca\n")
    for (k, v) in I.selected_sublevel_structure
        println(k, ": ", v)
    end
end

Base.show(io::IO, I::Ca40) = println(io, "⁴⁰Ca")  # suppress long output

Base.getindex(I::Ion, state::Union{Tuple{String,Real},String,Int}) = ionstate(I, state)

function Base.getproperty(I::Ion, s::Symbol)
    if s == :number || s == :position
        if typeof(getfield(I, s)) <: Missing
            @warn "ion has not been added to a configuration"
        return missing
        end
    end
    getfield(I, s)
end

# NEEDS TO BE CHANGED
function Base.setproperty!(I::Ion, s::Symbol, v::Tv) where{Tv}
    if (s == :mass || 
        s == :level_structure || 
        s == :shape || 
        s == :matrix_elements ||
        s == :selected_matrix_elements ||
        s == :number ||
        s == :position)
        return
    elseif s == :selected_level_structure
        @assert Tv == Vector{String} "type must be Vector{String}" 
        _, sls_dict, _, me_dict = _structure(v)
        Core.setproperty!(I, :selected_level_structure, sls_dict)
        Core.setproperty!(I, :selected_matrix_elements, me_dict)
        Core.setproperty!(I, :shape, [length(sls_dict)])
        I.stark_shift = OrderedDict{String,Real}()
        for key in v
            I.stark_shift[key] = 0.0
        end
        return
    end
    Core.setproperty!(I, s, v)
end

function Base.:(==)(b1::T, b2::T) where {T<:Ion}
    (
        b1.mass == b2.mass &&
        b1.selected_sublevel_structure == b2.selected_sublevel_structure &&
        b1.shape == b2.shape &&
        b1.stark_shift == b2.stark_shift
    )
end



#############################################################################################
# Add code for individual ion species
#############################################################################################

include("species/include_species.jl")