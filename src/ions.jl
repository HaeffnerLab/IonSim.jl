using WignerSymbols: wigner3j
using .PhysicalConstants: e, ca40_qubit_transition_frequency, m_ca40, ħ, α, μB


export mass, level_structure, selected_level_structure, stark_shift,
       selected_matrix_elements, matrix_element, get_basis, ion_number, ion_position,
       gJ, zeeman_shift, matrix_elements, zero_stark_shift, Ion, Ca40


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
level_structure(I::Ion)::OrderedDict{String,NamedTuple} = I.level_structure
selected_level_structure(I::Ion)::OrderedDict{String,NamedTuple} = I.selected_level_structure
matrix_elements(I::Ion)::OrderedDict = I.matrix_elements
selected_matrix_elements(I::Ion)::OrderedDict = I.selected_matrix_elements
ion_number(I::Ion)::Union{Int,Missing} = I.number
ion_position(I::Ion)::Union{Real,Missing} = I.position
stark_shift(I::Ion)::OrderedDict{String,Real} = I.stark_shift


#############################################################################################
# Ca40 Ion
#############################################################################################

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

function _structure(selected_level_structure)
    ED = ca40_qubit_transition_frequency
    fls = OrderedDict(
            "S-1/2" => (l=0, j=1//2, mⱼ=-1//2, E=0),
            "S+1/2" => (l=0, j=1//2, mⱼ=1//2, E=0),
            "D-5/2" => (l=2, j=5//2, mⱼ=-5//2, E=ED),
            "D-3/2" => (l=2, j=5//2, mⱼ=-3//2, E=ED),
            "D-1/2" => (l=2, j=5//2, mⱼ=-1//2, E=ED),
            "D+1/2" => (l=2, j=5//2, mⱼ=1//2, E=ED),
            "D+3/2" => (l=2, j=5//2, mⱼ=3//2, E=ED),
            "D+5/2" => (l=2, j=5//2, mⱼ=5//2, E=ED)
        )
    me = OrderedDict{Tuple,Function}()
    k = fls.keys
    for i in eachindex(k), j in i+1:length(k)
        t1, t2 = fls[k[i]], fls[k[j]]
        if ! (typeof(_ca40_matrix_elements((t1, t2), 0, 0, 0)) <: Nothing)
            f(Efield, γ, ϕ) = _ca40_matrix_elements((t1, t2), Efield, γ, ϕ)
            push!(me, (k[i], k[j]) => f)
        end
    end
    if selected_level_structure == "default" 
        sls = collect(keys(fls))
    else
        sls = selected_level_structure
    end
    sls_dict = OrderedDict{String,NamedTuple}()
    for k in sls
        @assert k in fls.keys "invalid level $k"
        push!(sls_dict, k => fls[k])
    end
    me_dict = OrderedDict{Tuple,Function}()
    for (k, v) in me
        if k[1] in sls && k[2] in sls
            push!(me_dict, k => v)
        end
    end
    return fls, sls_dict, me, me_dict
end

# geometric part of the matrix element for 40Ca S1/2 <-> D5/2 transitions, 
# assuming linearly polarized light
_ca40_geo = [
        (γ, ϕ) -> begin 
                γ = deg2rad(γ)
                ϕ = deg2rad(ϕ) 
                (1/2)abs(cos(γ)sin(2ϕ)) 
            end,
        (γ, ϕ) -> begin 
                γ = deg2rad(γ)
                ϕ = deg2rad(ϕ)
                sqrt(1/6)abs(cos(γ)cos(2ϕ) + im*sin(γ)cos(ϕ)) 
            end,
        (γ, ϕ) -> begin 
                γ = deg2rad(γ)
                ϕ = deg2rad(ϕ)
                sqrt(1/6)abs((1/2)cos(γ)sin(2ϕ) + im*sin(γ)sin(ϕ)) 
            end
    ]

function _matrix_elements(
        transition::Tuple{NamedTuple,NamedTuple}, Efield::Real, γ::Real, ϕ::Real
    )
    t1 = transition[1]
    t2 = transition[2]
    Δl = t2.l - t1.l
    Δj = t2.j - t1.j
    Δm = Int(abs(t2.m - t1.m))
    if abs(Δm) > 2
        return nothing
    end
    Ω_E1 = 0
    Ω_E2 = 0
    Ω_M1 = 0
    if (Δl ≡ 1 || Δl ≡ -1) && (abs(Δj) <= 1)
        #Ω_E1 = E1 matrix element
    end
    if (Δl ≡ 0 || Δl ≡ 2 || Δl ≡ -2) && (abs(Δj) <= 2)
        #Ω_E2 = E2 matrix element
    end
    if (Δl ≡ 0) && (abs(Δj) <= 1)
        #Ω_M1 = M1 matrix element
    end
    Ω_E1, Ω_E2, Ω_M1
end

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
function matrix_element(C::Ca40, transition::Vector{String}, Efield::Real, γ::Real, ϕ::Real)
    t1 = C.level_structure[transition[1]]
    t2 = C.level_structure[transition[2]]
    _ca40_matrix_elements((t1, t2), Efield, γ, ϕ)
end

function Base.print(I::Ca40)
    println("⁴⁰Ca\n")
    for (k, v) in I.selected_level_structure
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
        b1.selected_level_structure == b2.selected_level_structure &&
        b1.shape == b2.shape &&
        b1.stark_shift == b2.stark_shift
    )
end