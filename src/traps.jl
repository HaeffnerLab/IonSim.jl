using QuantumOptics: tensor, CompositeBasis
using .PhysicalConstants: ħ, c

export Trap, get_basis, Efield_from_pi_time, Efield_from_pi_time!, 
       transition_frequency, set_gradient!, Efield_from_rabi_frequency, 
       Efield_from_rabi_frequency!, global_beam!, get_η

   
#############################################################################################
# an ion trap with a linear ion chain configuration
#############################################################################################

"""
    Trap(;
            configuration::LinearChain, B::Real=0, Bhat::NamedTuple{(:x,:y,:z)=ẑ, ∇B::Real=0,
             δB::Union{Real,Function}=0, lasers::Vector{Laser}
    )

Information necessary to describe the Hamiltonian for a collection of ions in a linear chain
interacting with laser light.
#### user-defined fields
* `configuration<:LinearChain`
* `B`: A real value describing the mean magnitude of the B-field [Tesla].
* `Bhat::NamedTuple{(:x,:y,:z)}`: Describes the direction of the B-field (defaults to ẑ).
* `∇B`: Magnitude of the B-field gradient. We assume that the gradient always points along the
        z-direction. [Tesla / meter]
* `δB::Function`: Time-dependence of the B-field [Tesla]
* `lasers::Array{<:Laser}`: For each laser in the array, the pointing field should contain
        an array of `Tuple{Int,Real}`. The first element specifies the index of an ion
        in the `ions` field that the laser interacts with. The second element specifies a 
        scaling factor for the strength of that interaction (to be used, e.g., for 
        modeling cross-talk).
#### derived fields:
* `basis<:CompositeBasis`: The basis for describing the combined system, ions + vibrational
        modes. If constructing the Hamiltonian explictly (with [`hamiltonian`](@ref)), then
        the ordering of the basis is set, by convention, as 
        ``ion₁ ⊗ ion₂ ⊗ ... ⊗ ion_N ⊗ mode₁ ⊗ mode₂ ⊗ ... ⊗ mode_N``, where the ion bases are
        ordered according to the order in `T.configuration.ions` and the vibrational modes
        are ordered according to the order in 
        `[T.configuration.vibrational_modes.x, T.configuration.vibrational_modes.y, 
        T.configuration.vibrational_modes.z]`.       
    E.g. for:

    ```
    chain = LinearChain(ions=[C1, C2], com_frequencies=(x=2e6,y=2e6,z=1e6), 
    selected_modes=(x=[1, 2], y=[], z=[1]))
    ```

    The ordering of the basis would be

    `C1.basis ⊗ C2.basis ⊗ chain.vibrational_modes.x[1].basis 
    ⊗ chain.vibrational_modes.x[2].basis ⊗ chain.vibrational_modes.z[1].basis`
    
    Otherwise, the ordering is according to the form of the initial state used in the solver.
"""
mutable struct Trap
    configuration::LinearChain
    B::Real
    Bhat::NamedTuple{(:x,:y,:z)}
    ∇B::Real
    δB::Function
    lasers::Array{<:Laser}
    basis::CompositeBasis
    function Trap(;
            configuration::LinearChain, B=0, Bhat=ẑ, ∇B=0, δB=0, lasers=Laser[]
        )
        warn = nothing
        for i in 1:length(lasers)
            if length(lasers[i].pointing) == 0
                for n in eachindex(configuration.ions)
                    push!(lasers[i].pointing, (n, 1.0))
                end
            end
            for j in i+1:length(lasers)
                if lasers[j] ≡ lasers[i]
                    lasers[j] = copy(lasers[i])
                    if isnothing(warn)
                        warn = "Some lasers point to the same thing. Making copies."
                        @warn warn
                    end
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
            configuration.ions..., 
            configuration.vibrational_modes.x...,
            configuration.vibrational_modes.y...,
            configuration.vibrational_modes.z...,
        )
        new(configuration, B, Bhat, ∇B, δBt, lasers, basis) 
    end
end

function Base.setproperty!(T::Trap, s::Symbol, v)
    if s == :δB
        typeof(v) <: Number ? vt(t) = v : vt = v
        Core.setproperty!(T, s, vt)
    elseif s == :basis
        return
    elseif s == :configuration
        Core.setproperty!(T, s, v)
        basis = tensor(
            v.ions..., v.vibrational_modes.x..., v.vibrational_modes.y...,
            v.vibrational_modes.z...,
        )
        Core.setproperty!(T, :basis, basis)
    else
        Core.setproperty!(T, s, v)
    end
end

Base.show(io::IO, T::Trap) = print(io, "Trap")  # suppress long output


#############################################################################################
# general functions
#############################################################################################

"""
    get_basis(T::Trap)
Returns the composite basis describing the Hilbert space for `T`.
"""
function get_basis(T::Trap)::CompositeBasis 
    tensor(
            T.configuration.ions..., 
            T.configuration.vibrational_modes.x...,
            T.configuration.vibrational_modes.y...,
            T.configuration.vibrational_modes.z...,
        )
end 

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
Compute the E-field needed to get a certain `pi_time` with a certain `laser`-`ion` 
`transition`.

Alternatively, one may use
```
Efield_from_pi_time(
            pi_time::Real, T::Trap, laser_index::Int, ion_index::Int, 
            transition::Union{Tuple{String,String},Vector{<:String}}
        )
```
which is the same as 
`Efield_from_pi_time(pi_time, T.Bhat, T.lasers[laser_index], T.configuration.ions[ion_index],
transition)`
"""
function Efield_from_pi_time(
    pi_time::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion, 
    transition::Union{Tuple{String,String},Vector{<:String}}
)   
    p = laser.pointing
    (γ, ϕ) = map(x -> rad2deg(acos(ndot(Bhat, x))), [laser.ϵ, laser.k])
    s_indx = findall(x -> x[1] == ion.number, p)
    @assert length(s_indx) > 0 "This laser doesn't shine on this ion"
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
    Efield_from_pi_time!(
            pi_time::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion, 
            transition::Union{Tuple{String,String},Vector{<:String}}
        )
Same as [`Efield_from_pi_time`](@ref), but updates `laser[:E]` in-place.
"""
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

"""
    Efield_from_rabi_frequency(
        Ω::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion, 
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
Compute the Efield needed to get a certain rabi frequency `Ω` with a certain `laser`-`ion` 
`transition`.

Alternatively, one may use
```
Efield_from_rabi_frequency(
            Ω::Real, T::Trap, laser_index::Int, ion_index::Int, 
            transition::Union{Tuple{String,String},Vector{<:String}}
        )
```
which is the same as 
`Efield_from_rabi_frequency(pi_time, T.Bhat, T.lasers[laser_index], 
T.configuration.ions[ion_index], transition)`
"""
function Efield_from_rabi_frequency(
    Ω::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion, 
    transition::Union{Tuple{String,String},Vector{<:String}})
    Efield_from_pi_time(1 / 2Ω, Bhat, laser, ion, transition)
end

function Efield_from_rabi_frequency(
        Ω::Real, T::Trap, laser_index::Int, ion_index::Int, 
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
    Efield_from_rabi_frequency(
            Ω, T.Bhat, T.lasers[laser_index], T.configuration.ions[ion_index],
            transition
        )
end

"""
    Efield_from_rabi_frequency!(
        Ω::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion, 
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
Same as [`Efield_from_rabi_frequency`](@ref), but updates `laser[:E]` in-place.
"""
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
Compute the transition frequency of the `ion`'s selected transition under Bfield `B`.
Alternatively, one may use:
```
transition_frequency(
        T::Trap, ion::Ion, transition::Union{Tuple{String,String},Vector{<:String}}
    )
```
which is the same as `transition_frequency(T.B, ion, transition)` or
```
transition_frequency(
        T::Trap, ion_index::Int, transition::Union{Tuple{String,String},Vector{<:String}}
    )
```
which is the same as `transition_frequency(T.B, T.configuration.ions[ion_index], transition)`.
"""
function transition_frequency(
        B::Real, ion::Ion, transition::Union{Tuple{String,String},Vector{<:String}}
    )
    diff([map(x -> zeeman_shift(B, ion.selected_level_structure[x]), transition)...])[1]
end

function transition_frequency(
        T::Trap, ion::Ion, transition::Union{Tuple{String,String},Vector{<:String}}
    )
    B = T.B + T.∇B * ion.position
    transition_frequency(B, ion, transition)
end

function transition_frequency(
        T::Trap, ion_index::Int, transition::Union{Tuple{String,String},Vector{<:String}}
    )
    ion = T.configuration.ions[ion_index]
    B = T.B + T.∇B * ion.position
    transition_frequency(B, ion, transition)
end

"""
    set_gradient!(
            T::Trap, ion_indxs::Tuple{Int,Int}, transition::Tuple{String,String}, f::Real
        )
Sets the Bfield gradient in place to achieve a detuning `f` between the `transition` of two
ions, which are assumed to be of the same species. `ion_indxs` refer to the 
ordering of the ions in the chain.
"""
function set_gradient!(
        T::Trap, ion_indxs::Tuple{Int,Int}, transition::Tuple{String,String}, f::Real
    )
    ion1 = T.configuration.ions[ion_indxs[1]]
    ion2 = T.configuration.ions[ion_indxs[2]]
    separation = abs(ion1.position - ion2.position)
    E1, E2 = map(x -> zeeman_shift(1, ion1.selected_level_structure[x]), transition)
    T.∇B = f / (abs(E2 - E1) * separation)
end

# In QunatumOptics.jl, this method will return true whenever the shapes of b1 and b2 match,
# but we'd like to distinguish, i.e., between Ion ⊗ mode1 ⊗ mode2 and Ion ⊗ mode2 ⊗ mode1
# when mode1.N == mode2.N but mode1.axis ≠ mode2.axis.
function (QuantumOptics.:(==)(b1::T, b2::T) 
    where {T<:CompositeBasis{<:Vector{Int},<:Tuple{Vararg{<:IonSimBasis}}}})
N = length(b1.bases)
if N ≠ length(b2.bases)
    return false
end
for i in 1:N
    if !(b1.bases[i] == b2.bases[i])
        return false
    end
end
return true
end

"""
    get_η(V::VibrationalMode, L::Laser, I::Ion)
The Lamb-Dicke parameter: 
``|k|cos(\\theta)\\sqrt{\\frac{\\hbar}{2m\\nu}}`` 
for a given vibrational mode, ion and laser.
"""
function get_η(V::VibrationalMode, L::Laser, I::Ion; scaled=false)
    @fastmath begin
        k = 2π / L.λ
        scaled ? ν = 1 : ν = V.ν
        x0 = √(ħ / (2 * I.mass * 2π * ν))
        cosθ = ndot(L.k, V.axis)
        k * x0 * cosθ * V.mode_structure[I.number]
    end
end

