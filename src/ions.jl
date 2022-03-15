using WignerSymbols: wigner3j, wigner6j
using LinearAlgebra: cross
using .PhysicalConstants: e, ħ, α, μB, e, eye3, c_rank1, c_rank2


export Ion, speciesproperties, sublevels, sublevel_aliases, shape, stark_shift, ionnumber,
       ionposition, mass, charge, nuclearspin, zero_stark_shift!, set_stark_shift!,
       set_sublevel_alias!, ionlevels, quantumnumbers, landegf, zeeman_shift, energy,
       transitionfrequency, leveltransitions, subleveltransitions, einsteinA, lifetime,
       matrix_element




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
sublevels(I::Ion)::Vector{Tuple{String,Real}} = I.sublevels
sublevel_aliases(I::Ion)::Dict{String,Tuple} = I.sublevel_aliases
shape(I::Ion)::Vector{Int} = I.shape
stark_shift(I::Ion)::OrderedDict{Tuple,Real} = I.stark_shift
ionnumber(I::Ion)::Union{Int,Missing} = I.ionnumber
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

validatesublevel(I::Ion, sublevel::Tuple{String,Real}) = @assert sublevel in sublevels(I) "ion does not contain sublevel $sublevel"
validatesublevel(I::Ion, alias::String) = validatesublevel(I, alias2sublevel(I, alias))


""" 
This needs a docstring
"""
function zero_stark_shift!(I::Ion)
    for sublevel in keys(stark_shift(I))
        I.stark_shift[sublevel] = 0.0
    end
end


"""
This needs a docstring
"""
function set_stark_shift!(I::Ion, sublevel::Tuple{String,Real}, shift::Real)
    validatesublevel(I, sublevel)
    I.stark_shift[sublevel] = shift
end
function set_stark_shift!(I::Ion, stark_shift_dict::Dict)
    for sublevel in keys(stark_shift_dict)
        set_stark_shift!(I, sublevel, stark_shift_dict[sublevel])
    end
end


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
function set_sublevel_alias!(I::Ion, aliasdict::Dict{String,Tuple{String,Real}})
    for (alias, sublevel) in aliasdict
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
    levelstruct = speciesproperties(I).full_level_structure[sublevel[1]]
    names = (:n, :i, :s, :l, :j, :f, :m)
    values = (levelstruct.n, nuclearspin(I), 1//2, levelstruct.l, levelstruct.j, levelstruct.f, Rational(sublevel[2]))
    return NamedTuple{names}(values)
end
function quantumnumbers(I::Ion, level_or_alias::String)
    # If the second argument is a String, it could be either a level name or the alias of a sublevel
    if level_or_alias in ionlevels(I)
        # Second argument is a level name. Leave out the m quantum number
        levelstruct = speciesproperties(I).full_level_structure[level_or_alias]
        names = (:n, :i, :s, :l, :j, :f)
        values = (levelstruct.n, nuclearspin(I), 1//2, levelstruct.l, levelstruct.j, levelstruct.f)
        return NamedTuple{names}(values)
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
function stark_shift(I::Ion, sublevel::Tuple{String,Real})
    validatesublevel(I, sublevel)
    stark_shift(I)[sublevel]
end
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
function zeeman_shift(I::Ion, sublevel::Tuple{String,Real}, B::Real)
    validatesublevel(I, sublevel)
    properties = speciesproperties(I)
    if haskey(properties, :nonlinear_zeeman) && haskey(properties.nonlinear_zeeman, sublevel)
        nonlinear = properties.nonlinear_zeeman[sublevel](B)
    else
        nonlinear = 0.0
    end
    return zeeman_shift(B, landegf(I, sublevel[1]), sublevel[2]) + nonlinear
end
zeeman_shift(I::Ion, alias::String, B::Real) = zeeman_shift(I, alias2sublevel(I, alias), B)


"""
This needs a docstring
"""
function energy(I::Ion, sublevel::Tuple{String,Real}; B=0, ignore_starkshift=false)
    validatesublevel(I, sublevel)
    E0 = speciesproperties(I).full_level_structure[sublevel[1]].E
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
function leveltransitions(I::Ion)
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
function subleveltransitions(I::Ion)
    list = []
    for transition in leveltransitions(I)
        (L1, L2) = transition
        sublevels1 = [sublevel for sublevel in sublevels(I) if sublevel[1]==L1]
        sublevels2 = [sublevel for sublevel in sublevels(I) if sublevel[1]==L2]
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
function einsteinA(I::Ion, L1::String, L2::String)
    @assert (L1, L2) in leveltransitions(I) "invalid transition $L1 -> $L2"
    return speciesproperties(I).full_transitions[(L1, L2)].einsteinA
end
einsteinA(I::Ion, Lpair::Tuple) = einsteinA(I, Lpair[1], Lpair[2])


"""
This needs a docstring
"""
function transitionmultipole(I::Ion, L1::String, L2::String)
    @assert (L1, L2) in leveltransitions(I) "invalid transition $L1 -> $L2"
    return speciesproperties(I).full_transitions[(L1, L2)].multipole
end
transitionmultipole(I::Ion, Lpair::Tuple) = transitionmultipole(I, Lpair[1], Lpair[2])


"""
This needs a docstring
Computes lifetime based on species properties, *not* just the levels which are present in this ion
"""
function lifetime(I::Ion, level::String)
    @assert level in keys(speciesproperties(I).full_level_structure) "Ion species $(typeof(I)) does not contain level $level"
    totaltransitionrate = 0.0
    for (transition, info) in speciesproperties(I).full_transitions
        if transition[2] == level
            totaltransitionrate += info.einsteinA
        end
    end
    if totaltransitionrate == 0.0
        return Inf
    else
        return 1.0/totaltransitionrate
    end
end


"""
This needs a docstring
"""
function matrix_element(j1::Real, j2::Real, f1::Real, f2::Real, m1::Real, m2::Real, I::Real, ΔE::Real, A12::Real, multipole::String, Efield::Real, khat::NamedTuple, ϵhat::NamedTuple, Bhat::NamedTuple=(;z=1))
    # Level 1 *must* be the lower level and level 2 *must* be the upper level
    # Note that in this function, I is the nuclear spin, not an ion

    k = 2π*ΔE/c
    q = Int(m2-m1)
    
    Bhat_array = [Bhat.x, Bhat.y, Bhat.z]
    ϵhat_array = [ϵhat.x, ϵhat.y, ϵhat.z]
    khat_array = [khat.x, khat.y, khat.z]

    # Rotate unit vectors so that Bhat = ̂z
    if Bhat == ẑ
        R = eye3
    else
        a = cross(Bhat_array, [0, 0, 1])/norm(cross(Bhat_array, [0, 0, 1]))
        theta = acos(Bhat_array[2])
        amatrix = [0 -a[3] a[2]; a[3] 0 -a[1]; -a[2] a[1] 0]
        R = eye3 + sin(theta)*amatrix + (1-cos(theta))*amatrix^2    # Rotation matrix in axis-angle representation (axis=a, angle=theta)
    end
    ϵhat_rotated = R*ϵhat_array
    khat_rotated = R*khat_array

    if multipole=="E1"
        if abs(q) > 1
            return 0
        else
            hyperfine_factor = abs(sqrt((2*f1+1)*(2*f2+1)) * wigner6j(j2, I, f2, f1, 1, j1))
            geometric_factor = abs(sqrt(2j2+1) * wigner3j(f2, 1, f1, -m2, q, m1) * (transpose(c_rank1[q+2,:]) * ϵhat_rotated))
            units_factor = abs(e*Efield/(2ħ) * sqrt(3*A12/(α*c*k^3)))
        end
    elseif multipole=="E2"
        if abs(q) > 2
            return 0
        else
            hyperfine_factor = abs(sqrt((2*f1+1)*(2*f2+1)) * wigner6j(j2, I, f2, f1, 2, j1))
            geometric_factor = abs(sqrt(2j2+1) * wigner3j(f2, 2, f1, -m2, q, m1) * (transpose(khat_rotated) * c_rank2[:,:,q+3] * ϵhat_rotated))
            units_factor = abs(e*Efield/(2ħ) * sqrt(15*A12/(α*c*k^3)))
        end
    else
        @error "calculation of atomic transition matrix element for transition type $type not currently supported"
    end
    return units_factor * hyperfine_factor * geometric_factor / 2π
end
function matrix_element(I::Ion, transition::Tuple, Efield::Real, khat::NamedTuple, ϵhat::NamedTuple, Bhat::NamedTuple=(;z=1))
    E1 = energy(I, transition[1], ignore_starkshift=true)
    E2 = energy(I, transition[2], ignore_starkshift=true)
    @assert E2 > E1 "transition must be formatted (lower level, upper level)"
    qn1 = quantumnumbers(I, transition[1])
    qn2 = quantumnumbers(I, transition[2])
    A12 = einsteinA(I, transition[1][1], transition[2][1])
    multipole = transitionmultipole(I, transition[1][1], transition[2][1])
    matrix_element(qn1.j, qn2.j, qn1.f, qn2.f, qn1.m, qn2.m, nuclearspin(I), E2-E1, A12, multipole, Efield, khat, ϵhat, Bhat)
end



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
        @assert level ∉ [k[1] for k in keys(sublevels)] "multiple instances of level $level in ion constructor call"
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
            @assert (level, m) ∉ sublevels "repeated instance of sublevel $m in state $level"
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

# function Base.print(I::Ca40)
#     println("⁴⁰Ca\n")
#     for (k, v) in I.selected_sublevel_structure
#         println(k, ": ", v)
#     end
# end

# Base.show(io::IO, I::Ca40) = println(io, "⁴⁰Ca")  # suppress long output

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

function Base.setproperty!(I::Ion, s::Symbol, v::Tv) where{Tv}
    if (s == :species_properties || 
        s == :shape || 
        s == :number || 
        s == :position)
        return
    elseif s == :sublevels
        Core.setproperty!(I, :sublevels, _construct_sublevels(v, speciesproperties(I)))
        # Also update the stark shift dict as necessary; keep old stark shift values and assign zero stark shift to new sublevels
        starkshift_full_old = stark_shift(I)
        starkshift_full_new = OrderedDict{Tuple,Real}()
        for sublevel in sublevels(I)
            starkshift_full_new[sublevel] = (haskey(starkshift_full_old, sublevel) ? starkshift_full_old[sublevel] : 0.)
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

function Base.:(==)(b1::T, b2::T) where {T<:Ion}
    # Takes two ions to be equal if they are the same species, contain the same sublevels, and have the same stark shifts
    # Does not care about sublevel aliases or the ordering of sublevels
    (
        b1.species_properties == b2.species_properties &&
        sort(b1.sublevels) == sort(b2.sublevels) &&
        sort(b1.stark_shift) == sort(b2.stark_shift)
    )
end