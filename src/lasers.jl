using .PhysicalConstants: c, ϵ₀

export Laser,
    wavelength,
    intensity,
    detuning,
    polarization,
    wavevector,
    phase,
    pointing,
    wavelength!,
    intensity!,
    detuning!,
    polarization!,
    wavevector!,
    phase!,
    pointing!,
    efield

"""
    Laser(;λ=nothing, E=0, Δ=0, ϵ=(x̂+ŷ)/√2, k=ẑ, ϕ=0, pointing::Array{Tuple{Int,Real}})
        
The physical parameters defining laser light.
**args**
* `λ::Union{Real,Nothing}`: the wavelength of the laser in meters
* `I::Union{Function,Real}`: laser intensity in W/m²
* `Δ`: static detuning from f = c/λ in [Hz]
* `ϵ::NamedTuple`: (ϵ.x, ϵ.y, ϵ.z), polarization direction, requires norm of 1
* `k::NamedTuple`: (k.x, k.y, k.z), propagation direction, requires norm of 1
* `ϕ::Union{Function,Real}`: time-dependent phase. of course, this can also be used to model a 
    time-dependent detuning. Units are in radians. Note: if this is set to a function of time,
    then when constructing a Hamiltonian with the `hamiltonian` function, the units of time
    will be as specified by the `timescale` keyword argument.
* `pointing`: an array of `Tuple{Int,Real}` for describing ion-laser pointing configuration.
    (first element of the tuple is the index for an ion and the second element is the scaling
    factor for the laser's Efield which must be between 0 and 1).
"""
mutable struct Laser
    λ::Union{Real, Nothing}
    I::Function
    Δ::Real
    ϵ::NamedTuple{(:x, :y, :z)}
    k::NamedTuple{(:x, :y, :z)}
    ϕ::Function
    pointing::Vector
    function Laser(;
        λ = nothing,
        I::TI = 0,
        Δ = 0,
        ϵ = (x̂ + ŷ) / √2,
        k = ẑ,
        ϕ::Tϕ = 0,
        pointing = Array{Tuple{Int, <:Real}}(undef, 0)
    ) where {TI, Tϕ}
        rtol = 1e-6
        @assert isapprox(norm(ϵ), 1, rtol = rtol) "!(|ϵ| = 1)"
        @assert isapprox(norm(k), 1, rtol = rtol) "!(|k| = 1)"
        # @assert isapprox(ndot(ϵ, k), 0, rtol=rtol) "!(ϵ ⟂ k)"
        # Above commented out until we figure out a better place to put this warning
        a = pointing
        (ion_num, scaling) = map(x -> getfield.(a, x), fieldnames(eltype(a)))
        @assert length(ion_num) == length(unique(ion_num)) (
            "a laser is pointing at the same ion twice"
        )
        for s in scaling
            @assert 0 <= s <= 1 "must have s ∈ [0,1]"
        end
        TI <: Number ? It(t) = I : It = I
        Tϕ <: Number ? ϕt(t) = ϕ : ϕt = ϕ
        return new(λ, It, Δ, ϵ, k, ϕt, pointing)
    end
    # for copying
    Laser(λ, I, Δ, ϵ, k, ϕ, pointing) = new(λ, I, Δ, ϵ, k, ϕ, pointing)
end

#############################################################################################
# Object fields
#############################################################################################

wavelength(laser::Laser) = laser.λ
intensity(laser::Laser) = laser.I
detuning(laser::Laser) = laser.Δ
polarization(laser::Laser) = laser.ϵ
wavevector(laser::Laser) = laser.k
phase(laser::Laser) = laser.ϕ
pointing(laser::Laser) = laser.pointing


#############################################################################################
# Setters
#############################################################################################

function wavelength!(laser::Laser, λ::Real)
    laser.λ = λ
end

function intensity!(laser::Laser, I::Function)
    laser.I = I
end
function intensity!(laser::Laser, I::Real)
    laser.I = (t -> I)
end

function detuning(laser::Laser, Δ::Real)
    laser.Δ = Δ
end

function polarization!(laser::Laser, ϵ::NamedTuple{(:x, :y, :z)})
    rtol = 1e-6
    @assert isapprox(norm(ϵ), 1, rtol = rtol) "!(|ϵ| = 1)"
    # if ! isapprox(ndot(ϵ, laser.k), 0, rtol=rtol)
    #     @warn "!(ϵ ⟂ k)"
    # end 
    laser.ϵ = ϵ
end

function wavevector!(laser::Laser, k::NamedTuple{(:x, :y, :z)})
    rtol = 1e-6
    @assert isapprox(norm(k), 1, rtol = rtol) "!(|k| = 1)"
    # if ! isapprox(ndot(k, laser.ϵ), 0, rtol=rtol)
    #     @warn "!(ϵ ⟂ k)"
    # end
    laser.k = k
end

function phase!(laser::Laser, ϕ::Function)
    laser.ϕ = ϕ
end
function phase!(laser::Laser, ϕ::Real)
    laser.ϕ = (t -> ϕ)
end

function pointing!(laser::Laser, p::Vector{Tuple{T1, T2}} where T1<:Int where T2<:Real)
    (ion_num, scaling) = map(x -> getfield.(p, x), fieldnames(eltype(p)))
    @assert length(ion_num) == length(unique(ion_num)) (
        "a laser is pointing at the same ion twice"
    )
    for s in scaling
        @assert 0 <= s <= 1 "must have s ∈ [0,1]"
    end
    laser.pointing = p
end

efield(I::Real) = √(2I/(c*ϵ₀))
efield(laser::Laser) = t -> efield(intensity(laser)(t))

#############################################################################################
# Base functions
#############################################################################################

function Base.:(==)(L1::Laser, L2::Laser)
    for field in fieldnames(Laser)
        if getfield(L1, field) != getfield(L2, field)
            return false
        end
    end
    return true
end

function Base.print(L::Laser)
    println("λ: ", L.λ, " m")
    println("Δ: ", L.Δ, " Hz")
    println("ϵ̂: ", "(x=$(L.ϵ.x), y=$(L.ϵ.y), z=$(L.ϵ.z))")
    println("k̂: ", "(z=$(L.k.x), y=$(L.k.y), z=$(L.k.z))")
    println("I(t=0): ", "$(L.I(0.0)) W/m²")
    return println("ϕ(t=0): ", "$(L.ϕ(0.0)) ⋅ 2π") # Should the 2π be here?
end