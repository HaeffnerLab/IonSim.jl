using QuantumOptics: tensor, CompositeBasis


export Trap, trap, configuration, Bfield, Bhat, gradB, deltaB, get_lasers, get_basis
export Efield_from_pi_time, Efield_from_pi_time!, transition_frequency, set_gradient!
export Efield_from_rabi_frequency, Efield_from_rabi_frequency!, global_beam!


"""
    Trap
A container with all of the information necessary for simulation. (i.e. a collection of ions, 
modes, lasers and their physical relationships).
"""
abstract type Trap end

# required fields
configuration(T::Trap)::IonConfiguration = T.configuration
Bfield(T::Trap)::Union{Real,Function} = T.B
Bhat(T::Trap)::NamedTuple{(:x,:y,:z)} = T.Bhat
gradB(T::Trap)::Real = T.∇B
deltaB(T::Trap)::Function = T.δB
get_lasers(T::Trap)::Vector{Laser} = T.lasers
function get_basis(T::Trap)::CompositeBasis 
    tensor(
            [I.basis for I in T.configuration.ions]..., 
            T.configuration.vibrational_modes.x...,
            T.configuration.vibrational_modes.y...,
            T.configuration.vibrational_modes.z...,
        )
    # tensor(
    #         [I.basis for I in T.configuration.ions]..., 
    #         [FockBasis(V.N) for V in T.configuration.vibrational_modes.x]...,
    #         [FockBasis(V.N) for V in T.configuration.vibrational_modes.y]...,
    #         [FockBasis(V.N) for V in T.configuration.vibrational_modes.z]...,
    #     )
end 


#############################################################################################
# a trap with a linear ion chain configuration
#############################################################################################

"""
    trap(;
            configuration::linearchain, B::Real=0, Bhat::NamedTuple{(:x,:y,:z)=ẑ, ∇B::Real=0,
             δB::Union{Real,Function}=0, lasers::Vector{Laser}
    )

Information necessary to describe the Hamiltonian for a  collection of ions in a linear chain
interacting with laser light.
#### user-defined fields
* `configuration <: linearchain`
* `B`: either a real value describing the magnitude of the B-field or a function for
       describing its time dependence [Tesla]
* `Bhat::NamedTuple{(:x,:y,:z)}`: for descibing the direction of the B-field, defaults to ẑ.
* `∇B`: magnitude of the B-field gradient. We assume that the gradient always points along the
        z-direction. [Tesla / meter]
* `δB::Function`: time-dependence of the B-field [Tesla]
* `lasers::Array{<:Laser}`: For each laser in the array, the pointing field should contain
        an array of `Tuple{Int,Real}`. The first element specifies the index of an ion
        in the `ions` field that the laser interacts with. The second element specifies a 
        scaling factor for the strength of that interaction (to be used, e.g., for 
        modelling cross-talk).
#### derived fields:
* `basis::CompositeBasis`: The basis for describing the combined system, ions + vibrational
        modes. The ordering of the basis is set, by convention, as 
        ion₁ ⊗ ion₂ ⊗ ... ⊗ ionₙ ⊗ mode₁ ⊗ mode₂ ⊗ ... ⊗ modeₙ, where the ion bases are
        ordered according to the order in `T.configuration.ions` and the vibrational modes
        are ordered according to the order in 
        `[T.configuration.vibrational_modes.x, T.configuration.vibrational_modes.y, 
        T.configuration.vibrational_modes.z`].
        
    E.g. for:

    ```
    chain = linearchain(ions=[C1, C2], com_frequencies=(x=2e6,y=2e6,z=1e6), 
    selected_modes=(x=[1, 2], y=[], z=[1]))
    ```

    The ordering of the basis would be

    `C1.basis ⊗ C2.basis ⊗ chain.vibrational_modes.x[1].basis 
    ⊗ chain.vibrational_modes.x[2].basis ⊗ chain.vibrational_modes.z[1].basis`
"""
mutable struct trap <: Trap
    configuration::linearchain
    B::Real
    Bhat::NamedTuple{(:x,:y,:z)}
    ∇B::Real
    δB::Function
    lasers::Array{<:Laser}
    basis::CompositeBasis
    function trap(;
            configuration::linearchain, B=0, Bhat=(x=0, y=0, z=1), ∇B=0, δB=0, 
            lasers=Laser[]
        )
        for i in 1:length(lasers)
            if length(lasers[i].pointing) == 0
                for n in eachindex(configuration.ions)
                    push!(lasers[i].pointing, (n, 1.0))
                end
            end
            for j in 1:i-1
                if lasers[j] == lasers[i]
                    lasers[j] = copy(lasers[i])
                    @warn "Some lasers point to the same memory address. Making copies."
                    break
                end
            end
        end
        @assert isapprox(norm(Bhat), 1, rtol=1e-6) "!(|$Bhat| = 1)"
        for (li, l) in enumerate(lasers), p in l.pointing
            @assert p[1] <= length(configuration.ions) (
                """lasers[$li] points at configuration.ions[$(p[1])], but there are only
                 $(length(configuration.ions)) ions."""
            )
        end
        typeof(δB) <: Number ?  δBt(t) = δB : δBt = δB
        basis = tensor(
            [I.basis for I in configuration.ions]..., 
            configuration.vibrational_modes.x...,
            configuration.vibrational_modes.y...,
            configuration.vibrational_modes.z...,
        )
        new(configuration, B, Bhat, ∇B, δBt, lasers, basis) 
    end
end

function Base.setproperty!(T::trap, s::Symbol, v)
    if s == :δB
        typeof(v) <: Number ? vt(t) = v : vt = v
        Core.setproperty!(T, s, vt)
        return
    end
    Core.setproperty!(T, s, v)
end

# function Base.getindex(T::trap, S::String)
#     if typeof(T.configuration) == linearchain
#         v = T.configuration.vibrational_modes
#         V = vcat(T.configuration.ions, v.x, v.y, v.z)
#         for obj in V
#             if obj.label == S
#                 return obj
#             end
#         end
#     end
# end

Base.show(io::IO, T::trap) = print(io, "trap")  # suppress long output


#############################################################################################
# general functions
#############################################################################################

"""
    global_beam!(T::Trap, laser::Laser)
Set `laser` to shine with full intensity on all ions in `Trap`.
"""
function global_beam!(T::Trap, laser::Laser)
    for n in eachindex(T.configuration.ions)
        push!(laser.pointing, (n, 1.0))
    end
end

"""
    Efield_from_pi_time(
        pi_time::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion, 
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
Compute the Efield needed to get a certain pi_time with a certain laser-ion-transition.
"""
function Efield_from_pi_time(
    pi_time::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion, 
    transition::Union{Tuple{String,String},Vector{<:String}}
)
    (γ, ϕ) = map(x -> rad2deg(acos(ndot(Bhat, x))), [laser.ϵ, laser.k])
    s_indx = findall(x -> x[1] == ion.number, laser.pointing)
    s = laser.pointing[s_indx[1]][2]
    Ω = s * ion.selected_matrix_elements[tuple(transition...)](1, γ, ϕ)
    if Ω < 1e-15
        # even when coupling strength is zero, numerical error causes it to be finite
        # (on order 1e-16), this is a band-aid to prevent users from unknowingly setting 
        # the E-field to something absurd (like 1e20 V/m)
        return Inf
    end
    1 / (2Ω * pi_time) 
end

"""
    Efield_from_rabi_frequency(
        Ω::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion, 
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
Compute the Efield needed to get a certain rabi frequency with a certain laser-ion-transition.
"""
function Efield_from_rabi_frequency(
    Ω::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion, 
    transition::Union{Tuple{String,String},Vector{<:String}})
    Efield_from_pi_time(1 / 2Ω, Bhat, laser, ion, transition)
end

"""
    Efield_from_pi_time(
            pi_time::Real, T::Trap, laser_index::Int, ion_index::Int, 
            transition::Union{Tuple{String,String},Vector{<:String}}
        ) 
`laser_index` and `ion_index` refer to the position of the desired laser, ion in 
`T.lasers` and `T.configuration.ions`.
"""
function Efield_from_pi_time(
        pi_time::Real, T::Trap, laser_index::Int, ion_index::Int, 
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
    Efield_from_pi_time(
            pi_time, T.Bhat, T.lasers[laser_index], T.configuration.ions[ion_index],
            transition
        )
end

"""
    Efield_from_rabi_frequency(
            Ω::Real, T::Trap, laser_index::Int, ion_index::Int, 
            transition::Union{Tuple{String,String},Vector{<:String}}
        ) 
`laser_index` and `ion_index` refer to the position of the desired laser, ion in 
`T.lasers` and `T.configuration.ions`.
"""
function Efield_from_rabi_frequency(
        Ω::Real, T::Trap, laser_index::Int, ion_index::Int, 
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
    Efield_from_rabi_frequency(
            Ω, T.Bhat, T.lasers[laser_index], T.configuration.ions[ion_index],
            transition
        )
end

function Efield_from_pi_time!(
        pi_time::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion, 
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
    Efield::Float64 = Efield_from_pi_time(pi_time, Bhat, laser, ion, transition)    
    laser.E = t -> Efield
end

function Efield_from_pi_time!(
        pi_time::Real, T::Trap, laser_index::Int, ion_index::Int, 
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
    Efield::Float64 = Efield_from_pi_time(pi_time, T, laser_index, ion_index, transition)    
    T.lasers[laser_index].E = t -> Efield
end

function Efield_from_rabi_frequency!(
    Ω::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion, 
    transition::Union{Tuple{String,String},Vector{<:String}}
)
Efield::Float64 = Efield_from_rabi_frequency(Ω, Bhat, laser, ion, transition)    
laser.E = t -> Efield
end

function Efield_from_rabi_frequency!(
        Ω::Real, T::Trap, laser_index::Int, ion_index::Int, 
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
    Efield::Float64 = Efield_from_rabi_frequency(Ω, T, laser_index, ion_index, transition)    
    T.lasers[laser_index].E = t -> Efield
end

"""
    transition_frequency(
        B::Real, ion::Ion, transition::Union{Tuple{String,String},Vector{<:String}}
    )
Compute the transition frequency of the ion's selected transition under Bfield `B`.
"""
function transition_frequency(
        B::Real, ion::Ion, transition::Union{Tuple{String,String},Vector{<:String}}
    )
    diff([map(x -> zeeman_shift(B, ion.selected_level_structure[x]), transition)...])[1]
end

"""
    transition_frequency(
            T::Trap, ion::Ion, transition::Union{Tuple{String,String},Vector{<:String}}
        )
take B-field as given by `T.B`
"""
function transition_frequency(
        T::Trap, ion::Ion, transition::Union{Tuple{String,String},Vector{<:String}}
    )
    B = T.B + T.∇B * ion.position
    transition_frequency(B, ion, transition)
end

"""
    transition_frequency(
            T::Trap, ion_index::Int, transition::Union{Tuple{String,String},Vector{<:String}}
        )
take B-field as given by `T.B` and ion as given by `T.configuration.ions[ion_index]`.
"""
function transition_frequency(
        T::Trap, ion_index::Int, transition::Union{Tuple{String,String},Vector{<:String}}
    )
    ion = T.configuration.ions[ion_index]
    B = T.B + T.∇B * ion.position
    transition_frequency(B, ion, transition)
end

"""
    set_gradient!(
            T::trap, ion_indxs::Tuple{Int,Int}, transition::Tuple{String,String}, f::Real
        )
Sets the Bfield gradient in place to achieve a detuning `f` between the `transition` of two
ions, which are assumed to be of the same species. `ion_indxs` refer to the 
ordering of the ions in the chain.
"""
function set_gradient!(
        T::trap, ion_indxs::Tuple{Int,Int}, transition::Tuple{String,String}, f::Real
    )
    ion1 = T.configuration.ions[ion_indxs[1]]
    ion2 = T.configuration.ions[ion_indxs[2]]
    separation = abs(ion1.position - ion2.position)
    E1, E2 = map(x->zeeman_shift(1, ion1.selected_level_structure[x]), transition)
    T.∇B = f / (abs(E2 - E1) * separation)
end
