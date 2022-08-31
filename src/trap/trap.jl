using QuantumOptics: tensor, CompositeBasis
using .PhysicalConstants: ħ, c

export Trap, get_basis

"""
    Trap(;
            configuration::LinearChain, B::Real=0, Bhat::NamedTuple{(:x,:y,:z)=ẑ, ∇B::Real=0,
             δB::Union{Real,Function}=0, lasers::Vector{Laser}
    )

an ion trap with a linear ion chain configuration
Information necessary to describe the Hamiltonian for a collection of ions in a linear chain
interacting with laser light.
**user-defined fields**
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
**derived fields**
* `_cnst_δB::Bool`: A Boolean flag signifying whether or not `δB` is a constant function.
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
    Bhat::NamedTuple{(:x, :y, :z)}
    ∇B::Real
    δB::Function
    lasers::Array{<:Laser}
    basis::CompositeBasis
    _cnst_δB::Bool
    function Trap(;
        configuration::LinearChain,
        B = 0,
        Bhat = ẑ,
        ∇B = 0,
        δB::TδB = 0,
        lasers = Laser[]
    ) where {TδB}
        warn = nothing
        for i in 1:length(lasers)
            if length(lasers[i].pointing) == 0
                for n in eachindex(configuration.ions)
                    push!(lasers[i].pointing, (n, 1.0))
                end
            end
            for j in (i + 1):length(lasers)
                if lasers[j] ≡ lasers[i]
                    lasers[j] = copy(lasers[i])
                    if isnothing(warn)
                        warn = "Some lasers point to the same thing. Making copies."
                        @warn warn
                    end
                end
            end
        end
        @assert isapprox(norm(Bhat), 1, rtol = 1e-6) "!(|$Bhat| = 1)"
        for (li, l) in enumerate(lasers), p in l.pointing
            @assert p[1] <= length(configuration.ions) (
                """lasers[$li] points at configuration.ions[$(p[1])], but there are only
                 $(length(configuration.ions)) ions."""
            )
        end
        if TδB <: Number
            _cnst_δB = true
            δBt(t) = δB
        else
            _cnst_δB = false
            δBt = δB
        end
        basis = tensor(
            configuration.ions...,
            configuration.vibrational_modes.x...,
            configuration.vibrational_modes.y...,
            configuration.vibrational_modes.z...,
        )
        return new(configuration, B, Bhat, ∇B, δBt, lasers, basis, _cnst_δB)
    end
end

function Base.setproperty!(T::Trap, s::Symbol, v::Tv) where {Tv}
    if s == :δB
        if Tv <: Number
            _cnst_δB = true
            vt(t) = v
        else
            _cnst_δB = false
            vt = v
        end
        Core.setproperty!(T, s, vt)
        Core.setproperty!(T, :_cnst_δB, _cnst_δB)
    elseif s == :basis || s == :_cnst_δB
        return
    elseif s == :configuration
        Core.setproperty!(T, s, v)
        basis = tensor(
            v.ions...,
            v.vibrational_modes.x...,
            v.vibrational_modes.y...,
            v.vibrational_modes.z...,
        )
        Core.setproperty!(T, :basis, basis)
    else
        Core.setproperty!(T, s, v)
    end
end

Base.show(io::IO, T::Trap) = print(io, "Trap")  # suppress long output

# In QuantumOptics.jl, this method will return true whenever the shapes of b1 and b2 match,
# but we'd like to distinguish, i.e., between Ion ⊗ mode1 ⊗ mode2 and Ion ⊗ mode2 ⊗ mode1
# when mode1.N == mode2.N but mode1.axis ≠ mode2.axis.
function (
    QuantumOptics.:(==)(
        b1::T,
        b2::T
    ) where {T <: CompositeBasis{<:Vector{Int}, <:Tuple{Vararg{<:IonSimBasis}}}}
)
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
ionintrap(trap::Trap, ion::Ion)
Returns a boolean that indicates whether `ion` is actually in `trap`. Useful for checking if an error needs to be thrown.
"""
function ionintrap(trap::Trap, ion::Ion)
    return ion in ions(trap.configuration)
end

"""
get_basis(T::Trap)
Returns the composite basis describing the Hilbert space for `T`.
"""
function get_basis(T::Trap)::CompositeBasis
    return tensor(
        T.configuration.ions...,
        T.configuration.vibrational_modes.x...,
        T.configuration.vibrational_modes.y...,
        T.configuration.vibrational_modes.z...,
    )
end
