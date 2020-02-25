using QuantumOptics: FockBasis


export Vibration, vibrational_mode


"""
    Vibration
the physical parameters defining a vibrational mode of a linear chain of ions
"""
abstract type Vibration end

"""
    vibrational_mode(label::String, ν::Real, δν=0.0, mode_structure::Vector{Real}; N::Int=10)

#### user-defined fields
* `label`: convenience name
* `ν`: frequency [Hz]
* `mode_structure`: the normalized eigenvector describing the collective motion of the ions
        belonging to this mode
* `δν`: function describing time-dependent fluctuations of `ν`
* `N`: dimension of the Hilbert space
#### derived fields
* `basis<:QuantumOptics.FockBasis`: FockBasis representation of mode, cutoff``=N-1``.
* `axis`: the axis of vibration
"""
mutable struct vibrational_mode <: Vibration
    label::String
    ν::Real
    mode_structure::Vector{Real}
    δν::Function
    N::Int
    basis::QuantumOptics.FockBasis
    axis::NamedTuple{(:x,:y,:z)}
    function vibrational_mode(label, ν, mode_structure; δν=0.0, N=10, axis=ẑ)
        typeof(δν) <: Number ? δνt(t) = δν : δνt = δν
        new(label, ν, mode_structure, δνt, N, FockBasis(N-1), axis)
    end
end

# suppress long output
Base.show(io::IO, v::vibrational_mode) = print(io, "$(v.label)(ν=$(v.ν), axis=$(v.axis))")

function Base.setproperty!(V::Vibration, s::Symbol, v)
    if s == :mode_structure || s == :basis || s == :axis
        return
    end
    if s == :N
        @assert typeof(v) <: Int "N must be a positive integer"
        @assert v > 0 "N must be a positive integer"
        Core.setproperty!(V, :basis, FockBasis(v-1))
    elseif s == :δν
        typeof(v) <: Number ? vt(t) = v : vt = v
        Core.setproperty!(V, s, vt)
        return
    end
    Core.setproperty!(V, s, v)
end



