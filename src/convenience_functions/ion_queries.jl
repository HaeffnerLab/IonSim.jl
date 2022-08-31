export Efield_from_pi_time,
    Efield_from_pi_time!,
    Efield_from_rabi_frequency,
    Efield_from_rabi_frequency!,
    transitionfrequency,
    transitionwavelength,
    matrix_element,
    zeeman_shift

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

"""
    transitionfrequency(ion::Ion, transition::Tuple, T::Trap; ignore_starkshift=false)
Retuns The frequency of the transition `transition` including the Zeeman shift experienced by `ion` given its position in `T`.

One may alternatively replace `ion` with `ion_index::Int`, which instead specifies the index of the intended ion within `T`.
"""
transitionfrequency(ion::Ion, transition::Tuple, T::Trap; ignore_starkshift = false) =
    transitionfrequency(
        ion,
        transition;
        B = Bfield(T, ion),
        ignore_starkshift = ignore_starkshift
    )
transitionfrequency(ion_index::Int, transition::Tuple, T::Trap; ignore_starkshift = false) =
    transitionfrequency(
        T.configuration.ions[ion_index],
        transition,
        T;
        ignore_starkshift = ignore_starkshift
    )

"""
    transitionfrequency(ion::Ion, transition::Tuple, T::Trap; ignore_starkshift=false)
Retuns The wavelength of the transition `transition` including the Zeeman shift experienced by `ion` given its position in `T`.

One may alternatively replace `ion` with `ion_index`::Int, which instead specifies the index of the intended ion within `T`.
"""
transitionwavelength(ion::Ion, transition::Tuple, T::Trap; ignore_starkshift = false) =
    transitionwavelength(
        ion,
        transition;
        B = Bfield(T, ion),
        ignore_starkshift = ignore_starkshift
    )
transitionwavelength(
    ion_index::Int,
    transition::Tuple,
    T::Trap;
    ignore_starkshift = false
) = transitionwavelength(
    T.configuration.ions[ion_index],
    transition,
    T;
    ignore_starkshift = ignore_starkshift
)

#FIXME: This doc string is self-referential. 
"""
    matrix_element(I::Ion, transition::Tuple, T::Trap, laser::Laser, time::Real)
Calls `matrix_element(I::Ion, transition::Tuple, Efield::Real, khat::NamedTuple, ϵhat::NamedTuple, Bhat::NamedTuple=(;z=1))`
with `Efield`, `khat`, and `ϵhat` evaluated for `laser` at time `time`, and `Bhat` evaluated for `T`.

One may alternatively replace `ion` with `ion_index`::Int, which instead specifies the index of the intended ion within `T`.
"""
matrix_element(I::Ion, transition::Tuple, T::Trap, laser::Laser, time::Real) =
    matrix_element(I, transition, laser.E(time), laser.k, laser.ϵ, T.Bhat)
matrix_element(ion_index::Int, transition::Tuple, T::Trap, laser::Laser, time::Real) =
    matrix_element(T.configuration.ions[ion_index], transition, T, laser, time)

"""
    zeeman_shift(I::Ion, sublevel, T::Trap)
Calls `zeeman_shift(I::Ion, sublevel, B::Real)` with `B` evaluated for ion `I` in `T`.

One may alternatively replace `ion` with `ion_index`::Int, which instead specifies the index of the intended ion within `T`.
"""
zeeman_shift(I::Ion, sublevel::Union{Tuple{String, Real}, String}, T::Trap) =
    zeeman_shift(I, sublevel, Bfield(T, I))
zeeman_shift(ion_index::Int, sublevel::Union{Tuple{String, Real}, String}, T::Trap) =
    zeeman_shift(T.configuration.ions[ion_index], sublevel, T)
