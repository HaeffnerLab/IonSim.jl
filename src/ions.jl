using WignerSymbols: wigner3j
using .PhysicalConstants: e, ca40_qubit_transition_frequency, m_ca40, ħ, α, μB


export mass, full_level_structure, selected_sublevel_structure, sublevel_aliases,
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
full_level_structure(I::Ion)::Dict{String,NamedTuple} = I.full_level_structure
selected_sublevel_structure(I::Ion)::OrderedDict{Tuple,NamedTuple} = I.selected_sublevel_structure
sublevel_aliases(I::Ion)::Dict{String,Tuple} = I.sublevel_aliases
shape(I::Ion)::Vector{Int} = I.shape
full_transitions(I::Ion)::Dict{Tuple,Real} = I.full_transitions
selected_transitions(I::Ion)::Dict{Tuple,Real} = I.selected_transitions
stark_shift(I::Ion)::OrderedDict{Tuple,Real} = I.stark_shift
nonlinear_zeeman(I::Ion)::Dict{Tuple,Function} = I.nonlinear_zeeman
ion_number(I::Ion)::Union{Int,Missing} = I.number
ion_position(I::Ion)::Union{Real,Missing} = I.position

function _structure(selected_sublevels, full_level_structure, full_transitions)
    # First, construct the dictionary for selected_sublevel_structure
    selected_sublevel_structure = OrderedDict{Tuple{String,Rational},NamedTuple}()
    for manifold in selected_sublevels
        # Ensure that the string is a valid level
        level = manifold[1]
        @assert level in keys(full_level_structure) "invalid level $key"
        @assert level ∉ [k[1] for k in keys(selected_sublevel_structure)] "multiple instances of level $key in ion constructor call"
        level_structure = full_level_structure[level]

        # Add chosen sublevels
        sublevels = manifold[2]
        f = level_structure.f
        m_allowed = Array(-f:f)
        if sublevels == "all"
            sublevels = m_allowed
        elseif !(typeof(sublevels) <: Array)
            sublevels = [sublevels]
        end
        for m in sublevels
            m = Rational(m)
            @assert m in m_allowed "Zeeman sublevel m = $m not allowed for state $key with f = $f"
            @assert (key, m) ∉ keys(selected_sublevel_structure) "repeated instance of sublevel $m in state $key"
            selected_sublevel_structure[(key, m)] = (l=level_structure.l, j=level_structure.j, f=f, m=m, E=level_structure.E)
        end
    end

    # Then, construct the dictionary for selected_transitions
    selected_transitions = Dict{Tuple{String,String},Real}()
    levels = [manifold[1] for manifold in selected_sublevels]
    for (level_pair, value) in full_transitions
        if level_pair[1] in levels && level_pair[2] in levels
            selected_transitions[level_pair] = value
        end
    end

    return selected_sublevel_structure, selected_transitions
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
function alias2sublevel(I::Ion, alias)
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
(Don't forget that there's a method above for just stark_shift(I::Ion))
"""
stark_shift(I::Ion, sublevel::Tuple{String,Real}) = I.stark_shift[sublevel]
stark_shift(I::Ion, alias::String) = stark_shift(I, alias2sublevel(I, alias))

# geometric part of the matrix element for 40Ca S1/2 <-> D5/2 transitions, 
# assuming linearly polarized light
# _ca40_geo = [
#         (γ, ϕ) -> begin 
#                 γ = deg2rad(γ)
#                 ϕ = deg2rad(ϕ) 
#                 (1/2)abs(cos(γ)sin(2ϕ)) 
#             end,
#         (γ, ϕ) -> begin 
#                 γ = deg2rad(γ)
#                 ϕ = deg2rad(ϕ)
#                 sqrt(1/6)abs(cos(γ)cos(2ϕ) + im*sin(γ)cos(ϕ)) 
#             end,
#         (γ, ϕ) -> begin 
#                 γ = deg2rad(γ)
#                 ϕ = deg2rad(ϕ)
#                 sqrt(1/6)abs((1/2)cos(γ)sin(2ϕ) + im*sin(γ)sin(ϕ)) 
#             end
#     ]

# function _ca40_matrix_elements(
#         transition::Tuple{NamedTuple,NamedTuple}, Efield::Real, γ::Real, ϕ::Real
#     )
#     t1 = transition[1]
#     t2 = transition[2]
#     Δl = t2.l - t1.l
#     Δm = Int(abs(t2.mⱼ - t1.mⱼ))
#     if Δl ≡ 0 || abs(Δm) > 2 
#         return nothing
#     end
#     λ = c / ca40_qubit_transition_frequency
#     Ω = (e * Efield / ħ) * √(5λ^3 * (1/1.17) / (2 * π^3 * c * α)) / (2π)
#     wig = abs(wigner3j(t1.j, t2.j - t1.j, t2.j, -t1.mⱼ, t1.mⱼ - t2.mⱼ, t2.mⱼ))
#     Ω * _ca40_geo[Δm+1](γ, ϕ) * wig
# end

"""
    matrix_element(transition::Vector{String}, Efield::Real, γ::Real, ϕ::Real)

Computes the coupling strengths of the various S ⟷ D transitions in ⁴⁰Ca.
See e.g. page 30 of 
[Roos's thesis](https://quantumoptics.at/images/publications/dissertation/roos_diss.pdf).
Only considers linearly polarized light fields.

### args
* `C`: Ca40 ion
* `transition`: i.e. ["S-1/2", "D-1/2"]
* `Efield`: magnitude of the electric field at the position of the ion [V/m]
* `γ`: ``ϵ̂⋅B̂`` (angle between laser polarization and B-field) 
* `ϕ`: ``k̂⋅B̂`` (angle between laser k-vector and B-field)
"""
# function matrix_element(C::Ca40, transition::Vector{String}, Efield::Real, γ::Real, ϕ::Real)
#     t1 = C.level_structure[transition[1]]
#     t2 = C.level_structure[transition[2]]
#     _ca40_matrix_elements((t1, t2), Efield, γ, ϕ)
# end

"""
This needs a docstring
These are currently placeholder methods
"""
function matrix_element(Δl::Int, j1::Real, j2::Real, f1::Real, f2::Real, Δm::Int, ΔE::Real, A12::Real, Efield::Real, khat::NamedTuple, ϵhat::NamedTuple, Bhat::NamedTuple=(;ẑ=0))
    # Decide type of transition
    # Rotate unit vectors so that B is in z-direction?
    # Calculate matrix element
end
function matrix_element(I::Ion, transition, Efield::Real, khat::NamedTuple, ϵhat::NamedTuple, Bhat::NamedTuple=(;ẑ=0))
    # from I and transition determine Δl, j1, j2, f1, f2, Δm, ΔE, A12
    matrix_element(Δl, j1, j2, f1, f2, Δm, ΔE, A12, Efield, khat, ϵhat, Bhat)
end
function matrix_element(I::Ion, transition, laser::Laser, T::Trap)
    # Figure out Efield, khat, ϵhat from laser; Bhat from trap
    matrix_element(I, transition, Efield, khat, ϵhat, Bhat)
end



function Base.print(I::Ca40)
    println("⁴⁰Ca\n")
    for (k, v) in I.selected_sublevel_structure
        println(k, ": ", v)
    end
end

Base.show(io::IO, I::Ca40) = println(io, "⁴⁰Ca")  # suppress long output


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
gJ(l::Real, j::Real; s::Real=1/2) = 3/2 + (s * (s + 1) - l * (l + 1)) / (2j * (j + 1))

"""
    zeeman_shift(B::Real, l::Real, j::Real, mⱼ::Real)
``ΔE = (μ_B/ħ) ⋅ g_J(l, j) ⋅ B ⋅ mⱼ / 2π``
### args
* `B`: magnitude of B-field at ion
* `l`: orbital angular momentum quantum number
* `j`: total angular momentum quantum number
* `mⱼ`: projection of total angular momentum along quantization axis
"""
zeeman_shift(B::Real, l::Real, j::Real, mⱼ::Real) = zeeman_shift(B, (l=l, j=j, mⱼ=mⱼ))
zeeman_shift(B::Real, p::NamedTuple) = (μB/ħ) * gJ(p.l, p.j) * B * p.mⱼ / 2π
zeeman_shift(B::Real, p::Tuple) = zeeman_shift(B, (l=p[1], j=p[2], mⱼ=p[3]))
zeeman_shift(;B::Real, l::Real, j::Real, mⱼ::Real) = zeeman_shift(B, (l=l, j=j, mⱼ=mⱼ))

Base.getindex(I::Ion, state::String) = ionstate(I, state)
Base.getindex(I::Ion, state::Int) = ionstate(I, state)

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