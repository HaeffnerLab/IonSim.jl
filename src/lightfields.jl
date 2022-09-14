using .PhysicalConstants: c, ϵ₀

export Laser, Microwave

abstract type LightField end

function Base.:(==)(LF1::LightField, LF2::LightField)
    if typeof(LF1) != typeof(LF2)
        return false
    end
    for field in fieldnames(typeof(LF1))
        if getfield(LF1, field) != getfield(LF2, field)
            return false
        end
    end
    return true
end

function Base.print(LF::LightField)
    println("λ: ", LF.λ, " m")
    println("Δ: ", LF.Δ, " Hz")
    println("ϵ̂: ", "(x=$(LF.ϵ.x), y=$(LF.ϵ.y), z=$(LF.ϵ.z))")
    println("k̂: ", "(z=$(LF.k.x), y=$(LF.k.y), z=$(LF.k.z))")
    println("E(t=0): ", "$(LF.E(0.0)) V/m")
    return println("ϕ(t=0): ", "$(LF.ϕ(0.0)) ⋅ 2π")
end


"""
    Laser(;λ=nothing, E=0, Δ=0, ϵ=(x̂+ŷ)/√2, k=ẑ, ϕ=0, pointing::Array{Tuple{Int,Real}})
        
The physical parameters defining laser light.
**args**
* `λ::Union{Real,Nothing}`: the wavelength of the laser in meters
* `E::Union{Function,Real}`: magnitude of the E-field in V/m
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
mutable struct Laser <: LightField
    λ::Union{Real, Nothing}
    E::Function
    Δ::Real
    ϵ::NamedTuple{(:x, :y, :z)}
    k::NamedTuple{(:x, :y, :z)}
    ϕ::Function
    pointing::Vector
    function Laser(;
        λ = nothing,
        E::TE = 0,
        Δ = 0,
        ϵ = (x̂ + ŷ) / √2,
        k = ẑ,
        ϕ::Tϕ = 0,
        pointing = Array{Tuple{Int, <:Real}}(undef, 0)
    ) where {TE, Tϕ}
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
        TE <: Number ? Et(t) = E : Et = E
        Tϕ <: Number ? ϕt(t) = ϕ : ϕt = ϕ
        return new(λ, Et, Δ, ϵ, k, ϕt, pointing)
    end
    # for copying
    Laser(λ, E, Δ, ϵ, k, ϕ, pointing) = new(λ, E, Δ, ϵ, k, ϕ, pointing)
end

function Base.setproperty!(L::Laser, s::Symbol, v::Tv) where {Tv}
    rtol = 1e-6
    if s == :ϵ
        @assert isapprox(norm(v), 1, rtol = rtol) "!(|ϵ| = 1)"
        # if ! isapprox(ndot(L.k, v), 0, rtol=rtol)
        #     @warn "!(ϵ ⟂ k)"
        # end 
    elseif s == :k
        @assert isapprox(norm(v), 1, rtol = rtol) "!(|k| = 1)"
        # if ! isapprox(ndot(v, L.ϵ), 0, rtol=rtol)
        #     @warn "!(ϵ ⟂ k)"
        # end 
    elseif s == :pointing
        b = Tv <: Vector{Tuple{Int64, Float64}} || Tv <: Vector{Tuple{Int64, Int64}}
        @assert b "type != Vector{Tuple{Int,Real}}"
        (ion_num, scaling) = map(x -> getfield.(v, x), fieldnames(eltype(v)))
        @assert length(ion_num) == length(unique(ion_num)) (
            "a laser is pointing at the same ion twice"
        )
        for s in scaling
            @assert 0 <= s <= 1 "must have s ∈ [0,1]"
        end
    elseif s == :E || s == :ϕ
        Tv <: Number ? vt(t) = v : vt = v
        Core.setproperty!(L, s, vt)
        return
    end
    return Core.setproperty!(L, s, v)
end

# This could probably be made into one function with Unitfuls? 
laser_from_intensity(λ, i, Δ, ϵ, k, ϕ, pointing) = Laser(λ, sqrt(2i / (c*ϵ₀)), Δ, ϵ, k, ϕ, pointing)
laser_from_waist_power(λ, P, Δ, ϵ, k, ϕ, pointing, w₀) = laser_from_intensity(λ, 2P/(π*w₀^2), Δ, ϵ, k, ϕ, pointing)

"""
    Microwave(;λ=nothing, E=0, Δ=0, ϵ=(x̂+ŷ)/√2, k=ẑ, ϕ=0)
        
The physical parameters defining microwave radiation.
**args**
* `λ::Union{Real,Nothing}`: the wavelength of the microwave in meters
* `E::Union{Function,Real}`: magnitude of the E-field in V/m
* `Δ`: static detuning from f = c/λ in [Hz]
* `ϵ::NamedTuple`: (ϵ.x, ϵ.y, ϵ.z), polarization direction, requires norm of 1
* `k::NamedTuple`: (k.x, k.y, k.z), propagation direction, requires norm of 1
* `ϕ::Union{Function,Real}`: time-dependent phase. of course, this can also be used to model a 
    time-dependent detuning. Units are in radians. Note: if this is set to a function of time,
    then when constructing a Hamiltonian with the `hamiltonian` function, the units of time
    will be as specified by the `timescale` keyword argument.
"""
mutable struct Microwave <: LightField
    λ::Union{Real, Nothing}
    E::Function
    Δ::Real
    ϵ::NamedTuple{(:x, :y, :z)}
    k::NamedTuple{(:x, :y, :z)}
    ϕ::Function
    pointing::Array{Tuple{Int, <:Real}}(undef, 0)
    function Microwave(;
        λ = nothing,
        E::TE = 0,
        Δ = 0,
        ϵ = (x̂ + ŷ) / √2,
        k = ẑ,
        ϕ::Tϕ = 0,
    ) where {TE, Tϕ}
        rtol = 1e-6
        @assert isapprox(norm(ϵ), 1, rtol = rtol) "!(|ϵ| = 1)"
        @assert isapprox(norm(k), 1, rtol = rtol) "!(|k| = 1)"
        # @assert isapprox(ndot(ϵ, k), 0, rtol=rtol) "!(ϵ ⟂ k)"
        # Above commented out until we figure out a better place to put this warning
        TE <: Number ? Et(t) = E : Et = E
        Tϕ <: Number ? ϕt(t) = ϕ : ϕt = ϕ

        return new(λ, Et, Δ, ϵ, k, ϕt, [])
    end
    # for copying
    Microwave(λ, E, Δ, ϵ, k, ϕ) = new(λ, E, Δ, ϵ, k, ϕ)
end


function Base.setproperty!(MW::Microwave, s::Symbol, v::Tv) where {Tv}
    rtol = 1e-6
    if s == :ϵ
        @assert isapprox(norm(v), 1, rtol = rtol) "!(|ϵ| = 1)"
        # if ! isapprox(ndot(L.k, v), 0, rtol=rtol)
        #     @warn "!(ϵ ⟂ k)"
        # end 
    elseif s == :k
        @assert isapprox(norm(v), 1, rtol = rtol) "!(|k| = 1)"
        # if ! isapprox(ndot(v, L.ϵ), 0, rtol=rtol)
        #     @warn "!(ϵ ⟂ k)"
        # end 
    elseif s == :E || s == :ϕ
        Tv <: Number ? vt(t) = v : vt = v
        Core.setproperty!(MW, s, vt)
        return
    end
    return Core.setproperty!(MW, s, v)
end

# This could probably be made into one function with Unitfuls?
microwave_from_B(λ, B, Δ, ϵ, k, ϕ) = Microwave(λ, B*c, Δ, ϵ, k, ϕ)
microwave_from_intensity(λ, i, Δ, ϵ, k, ϕ, pointing) = Microwave(λ, sqrt(2i / (c*ϵ₀)), Δ, ϵ, k, ϕ)