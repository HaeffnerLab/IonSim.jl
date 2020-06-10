using QuantumOptics: Basis, basisstate
import Base: ==


export Vibration, vibrational_mode


"""
    Vibration
the physical parameters defining a vibrational mode of a linear chain of ions
"""
abstract type Vibration <: Basis end

"""
    vibrational_mode(
            ν::Real, mode_structure::Vector{Real}, δν::Union{Function,Real}=0.; N::Int=10,
            axis::NamedTuple{(:x,:y,:z)}=ẑ
        )

#### user-defined fields
* `ν::Real`: frequency [Hz]
* `mode_structure::Vector{Real}`: the normalized eigenvector describing the collective motion 
        of the ions belonging to this mode
* `δν::Union{Function,Real}`: either a function describing time-dependent fluctuations of `ν`
        or a real number which will be converted to the constant function `t -> δν`
* `N::Int`: highest modeled level of the corresponding basis
* `axis::NamedTuple{(:x,:y,:z)}`: the axis of symmetry for the vibration. this must lie along
        one of the basis vectors `x̂`, `ŷ` or `ẑ`
#### derived fields
* `shape::Vector{Int}`: indicates dimension of used Hilbert space (`=[N+1]`)

Note: the iᵗʰ Fock state (|i⟩) by indexing as `v=vibrational_mode(...); v[i]`
"""
mutable struct vibrational_mode <: Vibration
    ν::Real
    mode_structure::Vector{Real}
    δν::Function
    N::Int
    shape::Vector{Int}
    axis::NamedTuple{(:x,:y,:z)}
    function vibrational_mode(ν, mode_structure; δν=0., N=10, axis=ẑ)
        typeof(δν) <: Number ? δνt(t) = δν : δνt = δν
        new(ν, mode_structure, δνt, N, [N+1], axis)
    end
end

==(b1::T, b2::T) where {T<:vibrational_mode} = b1===b2

# suppress long output
Base.show(io::IO, v::vibrational_mode) = print(io, 
    "Vibration(ν=$(round(v.ν,sigdigits=4)), axis=$(_print_axis(v.axis)), N=$(v.N))")

function Base.setproperty!(V::Vibration, s::Symbol, v)
    if s == :mode_structure || s == :basis || s == :axis
        return
    end
    if s == :N
        @assert typeof(v) <: Int "N must be a positive integer"
        @assert v >= 0 "N must be a nonnegative integer"
        Core.setproperty!(V, :N, v)
        Core.setproperty!(V, :shape, Int[v+1])
    elseif s == :shape
        return
    elseif s == :δν
        typeof(v) <: Number ? vt(t) = v : vt = v
        Core.setproperty!(V, s, vt)
        return
    end
    Core.setproperty!(V, s, v)
end

function Base.getindex(V::vibrational_mode, n::Int)
    @assert 0 <= n <= V.N "n must be less than $(V.N+1)"
    basisstate(V, n+1)
end



