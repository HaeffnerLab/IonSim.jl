export matrix_element, lifetime, einsteinA, transitionmultipole
using LinearAlgebra: cross
using WignerSymbols: wigner3j, wigner6j

export energy,
    transitionfrequency,
    transitionmultipole,
    transitionwavelength,
    einsteinA,
    lifetime,
    matrix_element

# This function is written to be able to accept either a level or sublevel in the second argument
# Since both levels and aliases are strings, multidispatach can't tell the difference, so the second method distinguishes these cases with an if statement.
"""
    energy(I::Ion, sublevel; B=0, ignore_starkshift=false)
Returns energy of `sublevel` of `I`. A Zeeman shift may be included by setting the value of the magnetic field `B`. The Stark shift may be omitted by setting `ignore_starkshift=true`.
"""
function energy(I::Ion, sublevel::Tuple{String, Real}; B = 0, ignore_starkshift = false)
    validatesublevel(I, sublevel)
    E0 = speciesproperties(I).full_level_structure[sublevel[1]].E
    zeeman = zeeman_shift(I, sublevel, B)
    stark = (ignore_starkshift ? 0.0 : stark_shift(I, sublevel))
    return E0 + zeeman + stark
end

"""
    energy(I::Ion, level::trSing)
Returns the energy of `level` of `I`.
"""
function energy(I::Ion, level_or_alias::String; B = 0, ignore_starkshift = false)
    # If the second argument is a String, it could be either a level name or the alias of a sublevel
    if level_or_alias in levels(I)
        # Second argument is a level name. Return the bare energy of that level.
        return speciesproperties(I).full_level_structure[level_or_alias].E
    else
        # Second argument is a sublevel alias.
        return energy(
            I,
            alias2sublevel(I, level_or_alias),
            B = B,
            ignore_starkshift = ignore_starkshift
        )
    end
end

"""
    transitionfrequency(I::Ion, transition::Tuple; B=0, ignore_starkshift=false)
`transition` is a Tuple of two sublevels or levels.

Computes the absolute values of the difference in energies between `transition[1]` and `transition[2]`.

If between sublevels, then the Zeeman shift may be included by setting the value of the magnetic field `B`, and Stark shifts may be omitted by setting `ignore_starkshift=true`.
"""
function transitionfrequency(I::Ion, transition::Tuple; B = 0, ignore_starkshift = false)
    # Multidispatch of the function energy should make this work regardless of whether the transition is between levels or sublevels, and regardless of whether or not aliases are used
    return abs(
        energy(I, transition[1], B = B, ignore_starkshift = ignore_starkshift) -
        energy(I, transition[2], B = B, ignore_starkshift = ignore_starkshift)
    )
end

"""
    transitionwavelength(I::Ion, transition::Tuple; B=0, ignore_starkshift=false)
Returns the wavelength corresponding to `transitionfrequency(I::Ion, transition::Tuple; B=0, ignore_starkshift=false)`.
"""
function transitionwavelength(I::Ion, transition::Tuple; B = 0, ignore_starkshift = false)
    return c /
           transitionfrequency(I, transition, B = B, ignore_starkshift = ignore_starkshift)
end

"""
    einsteinA(I::Ion, Lpair::Tuple)
Returns Einstein A coefficient corresponding to the transition `Lpair[1] -> Lpair[2]`. The first level must be the lower level and the second must be the upper level.
"""
function einsteinA(I::Ion, L1::String, L2::String)
    @assert (L1, L2) in leveltransitions(I) "invalid transition $L1 -> $L2"
    return speciesproperties(I).full_transitions[(L1, L2)].einsteinA
end
einsteinA(I::Ion, Lpair::Tuple) = einsteinA(I, Lpair[1], Lpair[2])

"""
    transitionmultipole(I::Ion, Lpair::Tuple)
Returns the transition multiple (`'E1'`, `'E2'`, etc.) corresponding to the transition `Lpair[1] -> Lpair[2]`. The first level must be the lower level and the second must be the upper level.
"""
function transitionmultipole(I::Ion, L1::String, L2::String)
    @assert (L1, L2) in leveltransitions(I) "invalid transition $L1 -> $L2"
    return speciesproperties(I).full_transitions[(L1, L2)].multipole
end
transitionmultipole(I::Ion, Lpair::Tuple) = transitionmultipole(I, Lpair[1], Lpair[2])

"""
    lifetime(I::Ion, level::String)
Computes lifetime of `level` by summing the transition rates out of `level`.

The sum is taken over all levels that the species may have, rather than the levels present in the instance `I`.
"""
function lifetime(I::Ion, level::String)
    @assert level in keys(speciesproperties(I).full_level_structure) "Ion species $(typeof(I)) does not contain level $level"
    totaltransitionrate = 0.0
    for (transition, info) in speciesproperties(I).full_transitions
        if transition[2] == level
            totaltransitionrate += info.einsteinA
        end
    end
    if totaltransitionrate == 0.0
        return Inf
    else
        return 1.0 / totaltransitionrate
    end
end

"""
    matrix_element(I::Ion, transition::Tuple, Efield::Real, khat::NamedTuple, ϵhat::NamedTuple, Bhat::NamedTuple=(;z=1))
Computes the matrix elements (units of Hz) between two energy sublevels
**args**
* `I`: Ion undergoing transition
* `transition`: Tuple of sublevels (full names or aliases) between which the transition is being calculated. Must be formatted such that `energy(transition[2]) > energy(transition[1])`
* `Efield`: Amplitude of driving electric field
* `khat`: Unit vector of light wavevector
* `ϵhat`: Unit vector of light polarization
* `Bhat`: Unit vector of magnetic field
"""
function matrix_element(
    j1::Real,
    j2::Real,
    f1::Real,
    f2::Real,
    m1::Real,
    m2::Real,
    I::Real,
    ΔE::Real,
    A12::Real,
    multipole::String,
    Efield::Real,
    khat::NamedTuple,
    ϵhat::NamedTuple,
    Bhat::NamedTuple = (; z = 1)
)
    # Level 1 *must* be the lower level and level 2 *must* be the upper level
    # Note that in this function, I is the nuclear spin, not an ion

    k = 2π * ΔE / c
    q = Int(m2 - m1)

    Bhat_array = [Bhat.x, Bhat.y, Bhat.z]
    ϵhat_array = [ϵhat.x, ϵhat.y, ϵhat.z]
    khat_array = [khat.x, khat.y, khat.z]

    # Rotate unit vectors so that Bhat = ̂z
    if Bhat == ẑ
        R = eye3
    else
        a = cross(Bhat_array, [0, 0, 1]) / norm(cross(Bhat_array, [0, 0, 1]))
        theta = acos(Bhat_array[3])
        amatrix = [0 -a[3] a[2]; a[3] 0 -a[1]; -a[2] a[1] 0]
        R = eye3 + sin(theta) * amatrix + (1 - cos(theta)) * amatrix^2    # Rotation matrix in axis-angle representation (axis=a, angle=theta)
    end
    ϵhat_rotated = R * ϵhat_array
    khat_rotated = R * khat_array

    if multipole == "E1"
        if abs(q) > 1
            return 0
        else
            hyperfine_factor =
                abs(sqrt((2 * f1 + 1) * (2 * f2 + 1)) * wigner6j(j2, I, f2, f1, 1, j1))
            geometric_factor = abs(
                sqrt(2j2 + 1) *
                wigner3j(f2, 1, f1, -m2, q, m1) *
                (transpose(c_rank1[q + 2, :]) * ϵhat_rotated)
            )
            units_factor = abs(e * Efield / (2ħ) * sqrt(3 * A12 / (α * c * k^3)))
            return units_factor * hyperfine_factor * geometric_factor / 2π
        end
    elseif multipole == "E2"
        if abs(q) > 2
            return 0
        else
            hyperfine_factor =
                abs(sqrt((2 * f1 + 1) * (2 * f2 + 1)) * wigner6j(j2, I, f2, f1, 2, j1))
            geometric_factor = abs(
                sqrt(2j2 + 1) *
                wigner3j(f2, 2, f1, -m2, q, m1) *
                (transpose(khat_rotated) * c_rank2[:, :, q + 3] * ϵhat_rotated)
            )
            units_factor = abs(e * Efield / (2ħ) * sqrt(15 * A12 / (α * c * k^3)))
            return units_factor * hyperfine_factor * geometric_factor / 2π
        end
    else
        @error "calculation of atomic transition matrix element for transition type $multipole not currently supported"
    end
end
function matrix_element(
    I::Ion,
    transition::Tuple,
    Efield::Real,
    khat::NamedTuple,
    ϵhat::NamedTuple,
    Bhat::NamedTuple = (; z = 1)
)
    SL1 = transition[1]
    SL2 = transition[2]
    L1 = sublevel2level(I, SL1)
    L2 = sublevel2level(I, SL2)

    E1 = energy(I, L1)
    E2 = energy(I, L2)
    @assert E2 > E1 "transition must be formatted (lower level, upper level)"
    qn1 = quantumnumbers(I, SL1)
    qn2 = quantumnumbers(I, SL2)
    A12 = einsteinA(I, L1, L2)
    multipole = transitionmultipole(I, L1, L2)

    return matrix_element(
        qn1.j,
        qn2.j,
        qn1.f,
        qn2.f,
        qn1.m,
        qn2.m,
        nuclearspin(I),
        E2 - E1,
        A12,
        multipole,
        Efield,
        khat,
        ϵhat,
        Bhat
    )
end
