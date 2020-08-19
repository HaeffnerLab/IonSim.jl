using WignerSymbols: wigner3j
using .PhysicalConstants: e, ca40_qubit_transition_frequency, m_ca40, ħ, α, μB


export mass, nuclearspin, full_level_structure, selected_sublevel_structure, sublevel_aliases,
       full_transitions, selected_transitions, get_basis, stark_shift, ion_number, ion_position,
       set_sublevel_alias!, gJ, zeeman_shift, matrix_elements, zero_stark_shift, Ion


#############################################################################################
# Ion - the physical parameters defining an ion's structure
#############################################################################################

"""
    Ion
The physical parameters defining an isolated ion's internal structure.
"""
abstract type Ion <: IonSimBasis end

# required fields
mass(I::Ion)::Real = I.mass
charge(I::Ion)::Real = I.charge
nuclearspin(I::Ion)::Rational = I.nuclearspin
sublevel_structure(I::Ion)::OrderedDict{Tuple,NamedTuple} = I.sublevel_structure
sublevel_aliases(I::Ion)::Dict{String,Tuple} = I.sublevel_aliases
shape(I::Ion)::Vector{Int} = I.shape
transitions(I::Ion)::Dict{Tuple,Real} = I.transitions
stark_shift(I::Ion)::OrderedDict{Tuple,Real} = I.stark_shift
ion_number(I::Ion)::Union{Int,Missing} = I.number
ion_position(I::Ion)::Union{Real,Missing} = I.position
species_properties(I::Ion)::NamedTuple = I.species_properties

function _structure(selected_sublevels, properties)
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

    # Construct the dictionary for sublevel_structure
    sublevel_structure = OrderedDict{Tuple{String,Rational},NamedTuple}()
    for manifold in selected_sublevels
        # Ensure that the string is a valid level
        level = manifold[1]
        @assert level in keys(full_level_structure) "invalid level $level"
        @assert level ∉ [k[1] for k in keys(sublevel_structure)] "multiple instances of level $level in ion constructor call"
        level_structure = full_level_structure[level]

        # Add chosen sublevels
        sublevels = manifold[2]
        f = level_structure.f
        if haskey(properties, :gfactors) && haskey(properties.gfactors, level)
            gf = properties.gfactors[level]
        else
            gf = landegf(level_structure.l, level_structure.j, f, properties.nuclearspin)
        end
        m_allowed = Array(-f:f)
        if sublevels == "all"
            sublevels = m_allowed
        elseif !(typeof(sublevels) <: Array)
            sublevels = [sublevels]
        end
        for m in sublevels
            m = Rational(m)
            @assert m in m_allowed "Zeeman sublevel m = $m not allowed for state $level with f = $f"
            @assert (level, m) ∉ keys(sublevel_structure) "repeated instance of sublevel $m in state $level"
            # Here we add the key :nonlinear_zeeman to the sublevel structure if and only if a nonlinear zeeman function has been specified for this sublevel in the species
            if haskey(properties, :nonlinear_zeeman) && haskey(properties.nonlinear_zeeman, (level, m))
                sublevel_structure[(level, m)] = (l=level_structure.l, j=level_structure.j, f=f, m=m, E=level_structure.E, gf=gf, nonlinear_zeeman=properties.nonlinear_zeeman[(level, m)])
            else
                sublevel_structure[(level, m)] = (l=level_structure.l, j=level_structure.j, f=f, m=m, E=level_structure.E, gf=gf)
            end
        end
    end

    # Then, construct the dictionary for transitions
    transitions = Dict{Tuple{String,String},Real}()
    levels = [manifold[1] for manifold in selected_sublevels]
    for (level_pair, value) in properties.full_transitions
        if level_pair[1] in levels && level_pair[2] in levels
            transitions[level_pair] = value
        end
    end

    return sublevel_structure, transitions
end


"""
This needs a docstring
"""
function set_sublevel_alias!(I::Ion, sublevel::Tuple{String,Real}, alias::String)
    @assert sublevel in keys(I.selected_sublevel_structure) "ion does not have sublevel $sublevel"
    I.sublevel_aliases[alias] = sublevel
end
function set_sublevel_alias!(I::Ion, pairs::Vector{Tuple{Tuple{String,Real},String}})
    for (sublevel, alias) in pairs
        set_sublevel_alias!(I, sublevel, alias)
    end
end

"""
This needs a docstring
"""
function alias2sublevel(I::Ion, alias::String)
    all_aliases = I.sublevel_aliases
    @assert alias in keys(all_aliases) "no sublevel with alias $alias"
    return all_aliases[alias]
end

"""
This needs a docstring
"""
sublevel_structure(I::Ion, sublevel::Tuple{String,Real}) = I.selected_sublevel_structure[sublevel]
sublevel_structure(I::Ion, alias::String) = sublevel_structure(I, alias2sublevel(I, alias))

"""
This needs a docstring
"""
einsteinA(I::Ion, L1::Tuple{String,Real}, L2::Tuple{String,Real}) = I.selected_transitions[(L1, L2)]
einsteinA(I::Ion, L1::Tuple{String,Real}, L2::String) = I.selected_transitions[(L1, alias2sublevel(I, L2))]
einsteinA(I::Ion, L1::String, L2::Tuple{String,Real}) = I.selected_transitions[(alias2sublevel(I, L1), L2)]
einsteinA(I::Ion, L1::String, L2::String) = I.selected_transitions[(alias2sublevel(I, L1), alias2sublevel(I, L2))]

"""
This needs a docstring
(Don't forget that there's a method above for just stark_shift(I::Ion))
"""
stark_shift(I::Ion, sublevel::Tuple{String,Real}) = I.stark_shift[sublevel]
stark_shift(I::Ion, alias::String) = stark_shift(I, alias2sublevel(I, alias))

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
    sl1 = sublevel_structure(I, transition[1])
    sl2 = sublevel_structure(I, transition[2])
    matrix_element(sl2.l-sl2.l, sl1.j, sl2.j, sl1.f, sl2.f, sl2.m-sl1.m, abs(sl2.E-sl1.E), einsteinA(I, sl1, sl2), Efield, khat, ϵhat, Bhat)
end
matrix_element(I::Ion, transition::Tuple, T::Trap, laser::Laser) = matrix_element(I, transition, laser.E, laser.k, laser.ϵ, T.Bhat)



function Base.print(I::Ca40)
    println("⁴⁰Ca\n")
    for (k, v) in I.selected_sublevel_structure
        println(k, ": ", v)
    end
end

Base.show(io::IO, I::Ca40) = println(io, "⁴⁰Ca")  # suppress long output

Base.getindex(I::Ion, state::Union{Tuple{String,Real},String,Int}) = ionstate(I, state)

#############################################################################################
# general functions
#############################################################################################

"""
    gJ(l::real, j::real; s::Real=1/2)
Landé g-factor

### args
* `l`: orbital angular momentum quantum number
* `j`: total angular momentum quantum number
* `s`: spin angular momentum quantum number (defaults to 1/2)
"""
landegj(l::Real, j::Real, s::Real=1//2) = 3//2 + (s*(s+1) - l*(l+1)) / (2j*(j+1))
"""
This needs a docstring
"""
landegf(l::Real, j::Real, f::Real, i::Real, s::Real=1//2) = landegj(l, j, s)/2 * (1 + ((j*(j+1) - i*(i+1)) / (f*(f+1))))


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
function zeeman_shift(B::Real, I::Ion, sublevel::Union{Tuple{String,Real},String})
    structure = sublevel_structure(I, sublevel)
    nonlinear = (haskey(structure, :nonlinear_zeeman) ? structure.nonlinear_zeeman(B) : 0.)
    return zeeman_shift(B, structure.gf, structure.m) + nonlinear
end



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

function zero_stark_shift(I::Ion)
    for k in keys(I.stark_shift)
        I.stark_shift[k] = 0.0
    end
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

# Add code for individual ion species
include("species/include_species.jl")