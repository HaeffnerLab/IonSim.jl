export Efield_from_pi_time,
    Efield_from_pi_time!, Efield_from_rabi_frequency, Efield_from_rabi_frequency!

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
    pi_time::Real,
    Bhat::NamedTuple{(:x, :y, :z)},
    laser::Laser,
    ion::Ion,
    transition::Tuple
)
    p = laser.pointing
    s_indx = findall(x -> x[1] == ionnumber(ion), p)
    @assert length(s_indx) > 0 "This laser doesn't shine on this ion"
    s = p[s_indx[1]][2]
    Ω = s * matrix_element(ion, transition, 1, laser.k, laser.ϵ, Bhat)
    if Ω < 1e-15
        # even when coupling strength is zero, numerical error causes it to be finite
        # (on order 1e-16), this is a band-aid to prevent users from unknowingly setting
        # the E-field to something absurd (like 1e20 V/m)
        return Inf
    end
    return 1 / (2Ω * pi_time)
end

function Efield_from_pi_time(
    pi_time::Real,
    T::Trap,
    laser_index::Int,
    ion_index::Int,
    transition::Tuple
)
    return Efield_from_pi_time(
        pi_time,
        T.Bhat,
        T.lasers[laser_index],
        T.configuration.ions[ion_index],
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
    pi_time::Real,
    Bhat::NamedTuple{(:x, :y, :z)},
    laser::Laser,
    ion::Ion,
    transition::Tuple
)
    Efield::Float64 = Efield_from_pi_time(pi_time, Bhat, laser, ion, transition)
    return laser.E = t -> Efield
end

function Efield_from_pi_time!(
    pi_time::Real,
    T::Trap,
    laser_index::Int,
    ion_index::Int,
    transition::Tuple
)
    Efield::Float64 = Efield_from_pi_time(pi_time, T, laser_index, ion_index, transition)
    return T.lasers[laser_index].E = t -> Efield
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
    Ω::Real,
    Bhat::NamedTuple{(:x, :y, :z)},
    laser::Laser,
    ion::Ion,
    transition::Tuple
)
    return Efield_from_pi_time(1 / 2Ω, Bhat, laser, ion, transition)
end

function Efield_from_rabi_frequency(
    Ω::Real,
    T::Trap,
    laser_index::Int,
    ion_index::Int,
    transition::Tuple
)
    return Efield_from_rabi_frequency(
        Ω,
        T.Bhat,
        T.lasers[laser_index],
        T.configuration.ions[ion_index],
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
    Ω::Real,
    Bhat::NamedTuple{(:x, :y, :z)},
    laser::Laser,
    ion::Ion,
    transition::Tuple
)
    Efield::Float64 = Efield_from_rabi_frequency(Ω, Bhat, laser, ion, transition)
    return laser.E = t -> Efield
end

function Efield_from_rabi_frequency!(
    Ω::Real,
    T::Trap,
    laser_index::Int,
    ion_index::Int,
    transition::Tuple
)
    Efield::Float64 = Efield_from_rabi_frequency(Ω, T, laser_index, ion_index, transition)
    return T.lasers[laser_index].E = t -> Efield
end
