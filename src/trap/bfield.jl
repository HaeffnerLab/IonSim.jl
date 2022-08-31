export transition_frequency

"""
Bfield(T::Trap, ion::Ion)
Retuns the value of the magnetic field in `T` at the location of `ion`, including both the trap's overall B-field and its B-field gradient.
"""
function Bfield(T::Trap, ion::Ion)
    @assert ionintrap(T, ion) "trap does not contain ion"
    return T.B + T.∇B * ionposition(ion)
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

"""
set_gradient!(
        T::Trap, ion_indxs::Tuple{Int,Int}, transition::Tuple, f::Real
    )
Sets the Bfield gradient in place to achieve a detuning `f` between the `transition` of two
ions, which are assumed to be of the same species. `ion_indxs` refer to the
ordering of the ions in the chain.
"""
function set_gradient!(T::Trap, ion_indxs::Tuple{Int, Int}, transition::Tuple, f::Real)
    ionA = T.configuration.ions[ion_indxs[1]]
    ionB = T.configuration.ions[ion_indxs[2]]
    separation = abs(ionposition(ionA) - ionposition(ionB))

    (SL1, SL2) = transition
    L1 = sublevel2level(ionA, SL1)
    L2 = sublevel2level(ionA, SL2)
    g1 = landegf(ionA, L1)
    g2 = landegf(ionA, L2)
    m1 = quantumnumbers(ionA, SL1).m
    m2 = quantumnumbers(ionA, SL2).m
    # Calculate Zeeman shifts with a unit B-field using a method of zeeman_shift that ensures a nonlinear term is not used
    E1 = zeeman_shift(1.0, g1, m1)
    E2 = zeeman_shift(1.0, g2, m2)
    return T.∇B = f / (abs(E2 - E1) * separation)
end
