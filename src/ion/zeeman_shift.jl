using IonSim.PhysicalConstants

export landegj, landegf, zeeman_shift
"""
    landegj(l::Real, j::Real, s::Real=1//2)
Landé g-factor of fine structure energy level

**args**
* `l`: orbital angular momentum quantum number
* `j`: electron total angular momentum quantum number
* `s`: electronic spin angular momentum quantum number (defaults to 1/2)
"""
landegj(l::Real, j::Real, s::Real = 1 // 2) =
    3 // 2 + (s * (s + 1) - l * (l + 1)) / (2j * (j + 1))

"""
    landegf(l::Real, j::Real, f::Real, i::Real, s::Real=1//2)
Landé g-factor of hyperfine energy level

**args**
* `l`: orbital angular momentum quantum number
* `j`: electron total angular momentum quantum number
* `f`: total angular momentum quantum number
* `i`: nuclear spin angular momentum quantum number
* `s`: electronic spin angular momentum quantum number (defaults to 1/2)
"""
landegf(l::Real, j::Real, f::Real, i::Real, s::Real = 1 // 2) =
    landegj(l, j, s) / 2 * (1 + ((j * (j + 1) - i * (i + 1)) / (f * (f + 1))))
landegf(qnums::NamedTuple) = landegf(qnums.l, qnums.j, qnums.f, qnums.i, qnums.s)

"""
    landegf(I::Ion, level::String)
`landegf` for the quantum numbers of `level` in `I`.
"""
function landegf(I::Ion, level::String)
    properties = speciesproperties(I)
    if !ismissing(properties.gfactors) && haskey(properties.gfactors, level)
        return properties.gfactors[level]
    else
        return landegf(quantumnumbers(I, level))
    end
end

"""
    zeeman_shift(I::Ion, sublevel}, B::Real)
Returns the Zeeman shift at a magnetic field of `B` of `sublevel` of `I`.

If `sublevel` has a custom g-factor defined, then this is used. Otherwise, `landegf` is used to compute the Landé g-factor.

Zeeman shift calculated as ``ΔE = (μ_B/ħ) ⋅ g_f ⋅ B ⋅ m / 2π``
"""
function zeeman_shift(I::Ion, sublevel::Tuple{String, Real}, B::Real)
    validatesublevel(I, sublevel)
    properties = speciesproperties(I)
    if !ismissing(properties.nonlinear_zeeman) &&
       haskey(properties.nonlinear_zeeman, sublevel)
        nonlinear = properties.nonlinear_zeeman[sublevel](B)
    else
        nonlinear = 0.0
    end
    return zeeman_shift(B, landegf(I, sublevel[1]), sublevel[2]) + nonlinear
end
zeeman_shift(B::Real, g::Real, m::Real) = (μB / ħ) * g * B * m / 2π
zeeman_shift(B::Real, l::Real, j::Real, f::Real, m::Real, i::Real, s::Real = 1 // 2) =
    zeeman_shift(B, landegf(l, j, f, i, s), m)
zeeman_shift(B::Real, qnums::NamedTuple) =
    zeeman_shift(B, qnums.l, qnums.j, qnums.f, qnums.m, qnums.i, qnums.s)
zeeman_shift(I::Ion, alias::String, B::Real) = zeeman_shift(I, alias2sublevel(I, alias), B)
