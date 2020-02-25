export Laser, laser

"""
    Laser
the physical parameters defining laser light
"""
abstract type Laser end

"""
    laser(;
            label="", 
            E=1e5, Δ=0, ϵ=(x=1/√2, y=1/√2, z=0), k=(x=0, y=0, z=1), ϕ=0, λ=729.147e-9,
            pointing::Array{Tuple{Int,Real}}
        )

### args
* `label`: convenience label
* `E::Union{Function,Real}`: magnitude of the E-field in V/m
* `Δ`: static detuning from f = c/λ in [Hz]
* `ϵ::NamedTuple`: (ϵ.x, ϵ.y, ϵ.z), polarization direction, requires norm of 1
* `k::NamedTuple`: (k.x, k.y, k.z), propagation direction, requires norm of 1
* `ϕ::Union{Function,Real}`: time-dependent phase. of course, this can also be used to model a 
    time-dependent detuning. units are Hz/s or multiples of π (depending on your perspective).
* `λ::Real`: the wavelength of the laser in meters
* `pointing`: an array of `Tuple{Int,Real}` for describing ion-laser pointing configuration.
    (first element of the tuple is the index for an ion and the second element is the scaling
    factor for the laser's Efield which must be between 0 and 1).
"""
mutable struct laser <: Laser
    label::String
    E::Function
    Δ::Real
    ϵ::NamedTuple{(:x,:y,:z)}
    k::NamedTuple{(:x,:y,:z)}
    ϕ::Function 
    λ::Real
    pointing::Vector
    function laser(; 
            label="", E=1e5, Δ=0, ϵ=(x=1/√2, y=1/√2, z=0), k=(x=0, y=0, z=1), ϕ=0, 
            λ=c/ca40_qubit_transition_frequency, pointing=Array{Tuple{Int,<:Real}}(undef, 0)
        )
        rtol = 1e-6
        @assert isapprox(norm(ϵ), 1, rtol=rtol) "!(|ϵ| = 1)"
        @assert isapprox(norm(k), 1, rtol=rtol) "!(|k| = 1)"
        @assert isapprox(ndot(ϵ, k), 0, rtol=rtol) "!(ϵ ⟂ k)"
        a = pointing
        (ion_num, scaling) = map(x->getfield.(a, x), fieldnames(eltype(a)))
        @assert length(ion_num) == length(unique(ion_num)) (
                "a laser is pointing at the same ion twice"
        )
        for s in scaling
            @assert 0 <= s <= 1 "must have s ∈ [0,1]"
        end
        typeof(E) <: Number ?  Et(t) = E : Et = E
        typeof(ϕ) <: Number ? ϕt(t) = ϕ : ϕt = ϕ
        new(label, Et, Δ, ϵ, k, ϕt, λ, pointing)
    end
    # for copying
    laser(label, E, Δ, ϵ, k, ϕ, λ, pointing) = new(label, E, Δ, ϵ, k, ϕ, λ, pointing)
end 

function Base.print(L::Laser)
    print("λ: ", L.λ, " nm\n")
    print("Δ: ", L.Δ, " Hz\n")
    print("ϵ̂: ", "(x=$(L.ϵ.x), y=$(L.ϵ.y), z=$(L.ϵ.z))\n")
    print("k̂: ", "(z=$(L.k.x), y=$(L.k.y), z=$(L.k.z))\n")
    print("E(t=0): ", "$(L.E(0.0)) V/m\n")
    print("ϕ(t=0): ", "$(L.ϕ(0.0)) ⋅ π\n")
end

function Base.setproperty!(L::laser, s::Symbol, v)
    rtol = 1e-6
    if s == :ϵ
        @assert isapprox(norm(v), 1, rtol=rtol) "!(|ϵ| = 1)"
        if ! isapprox(ndot(L.k, v), 0, rtol=rtol)
            @warn "!(ϵ ⟂ k)"
        end 
    elseif s == :k
        @assert isapprox(norm(v), 1, rtol=rtol) "!(|k| = 1)"
        if ! isapprox(ndot(v, L.ϵ), 0, rtol=rtol)
            @warn "!(ϵ ⟂ k)"
        end 
    elseif s == :pointing
        b = typeof(v)<:Vector{Tuple{Int64,Float64}} || typeof(v)<:Vector{Tuple{Int64,Int64}}
        @assert b  "type != Vector{Tuple{Int,Real}}"
        (ion_num, scaling) = map(x->getfield.(v, x), fieldnames(eltype(v)))
        @assert length(ion_num) == length(unique(ion_num)) (
                "a laser is pointing at the same ion twice"
        )
        for s in scaling
            @assert 0 <= s <= 1 "must have s ∈ [0,1]"
        end
    elseif s == :E || s == :ϕ 
        typeof(v) <: Number ? vt(t) = v : vt = v
        Core.setproperty!(L, s, vt)
        return
    end
    Core.setproperty!(L, s, v)
end
