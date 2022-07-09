using WignerSymbols: wigner3j, wigner6j
using LinearAlgebra: cross
#using .PhysicalConstants: e, ħ, α, μB, e, eye3, c_rank1, c_rank2
using IonSim.PhysicalConstants

export Ion,
    IonProperties,
    speciesproperties,
    sublevels,
    sublevel_aliases,
    sublevelalias,
    shape,
    stark_shift,
    ionnumber,
    ionposition,
    mass,
    charge,
    nuclearspin,
    zero_stark_shift!,
    set_stark_shift!,
    alias2sublevel,
    sublevel2level,
    set_sublevel_alias!,
    clear_sublevel_alias!,
    clear_all_sublevel_aliases!,
    levels,
    quantumnumbers,
    landegf,
    zeeman_shift,
    energy,
    transitionfrequency,
    transitionwavelength,
    leveltransitions,
    subleveltransitions,
    einsteinA,
    transitionmultipole,
    lifetime,
    matrix_element,
    IonInstance

#############################################################################################
# Ion - the physical parameters defining an ion's structure
#############################################################################################

"""
    Ion
The physical parameters defining an isolated ion's internal structure.
"""
abstract type Ion <: IonSimBasis end

#############################################################################################
# Object fields
#############################################################################################

speciesproperties(I::Ion)::IonProperties = I.species_properties
sublevels(I::Ion)::Vector{Tuple{String, Real}} = I.sublevels
sublevel_aliases(I::Ion)::Dict{String, Tuple} = I.sublevel_aliases
shape(I::Ion)::Vector{Int} = I.shape
stark_shift(I::Ion)::OrderedDict{Tuple, Real} = I.stark_shift
ionnumber(I::Ion)::Union{Int, Missing} = I.ionnumber
ionposition(I::Ion)::Union{Real, Missing} = I.position

#############################################################################################
# General properties of species
#############################################################################################

mass(x) = print(x)
mass(I::Ion)::Real = speciesproperties(I).mass
charge(I::Ion)::Real = speciesproperties(I).charge * e
nuclearspin(I::Ion)::Rational = speciesproperties(I).nuclearspin

#############################################################################################
# Functions to modify ion properties
#############################################################################################

validatesublevel(I::Ion, sublevel::Tuple{String, Real}) =
    @assert sublevel in sublevels(I) "Ion does not contain sublevel $sublevel. Use sublevels(Ion) to see list of available sublevels."
validatesublevel(I::Ion, alias::String) = validatesublevel(I, alias2sublevel(I, alias))

"""
    set_sublevel_alias!(I::Ion, sublevel::Tuple{String,Real}, alias::String)
Assigns an alias `alias` to `sublevel` of `I`. This then allows one to pass `alias` in place of `sublevel` (for `I` only) into any function which accepts a sublevel as an argument.
"""
function set_sublevel_alias!(I::Ion, sublevel::Tuple{String, Real}, alias::String)
    validatesublevel(I, sublevel)
    @assert alias ∉ levels(I) "cannot make alias name identical to level name ($alias)"
    sublevel_rational = (sublevel[1], Rational(sublevel[2]))   # Force m to be Rational
    return I.sublevel_aliases[alias] = sublevel_rational
end
"""
    set_sublevel_alias!(I::Ion, aliasassignments)
`aliasassignments` specifies which aliases should be assigned to which sublevels. There are two options to do this:
* `aliasassignments` is a Vector of Tuples, with the first element of each being the sublevel (`sublevel::Tuple{String,Real}`) and the second being its assigned alias (`alias::String`)
* `aliasassignments` is a Dict with the format `alias::String => sublevel::Tuple{String,Real}`
Calls `set_sublevel_alias!(I, sublevel, alias)` for each pair `sublevel, alias`.
"""
function set_sublevel_alias!(
    I::Ion,
    pairs::Vector{Tuple{Tuple{String, R}, String}} where {R <: Real}
)
    for (sublevel, alias) in pairs
        set_sublevel_alias!(I, sublevel, alias)
    end
end
function set_sublevel_alias!(
    I::Ion,
    aliasdict::Dict{String, Tuple{String, R}} where {R <: Real}
)
    for (alias, sublevel) in aliasdict
        set_sublevel_alias!(I, sublevel, alias)
    end
end

"""
    clear_sublevel_alias!(I::Ion, sublevel)
Erases the assignment of an alias to `sublevel` of Ion `I`. Accepts either the full sublevel `Tuple{String,Real}` or its alias `String`.
Also accepts a vector of sublevels to clear multiple alias assignments in a single call.
"""
function clear_sublevel_alias!(I::Ion, sublevel::Tuple{String, Real})
    alias = sublevelalias(I, sublevel)
    return delete!(I.sublevel_aliases, alias)
end
function clear_sublevel_alias!(I::Ion, alias::String)
    return delete!(I.sublevel_aliases, alias)
end
function clear_sublevel_alias!(I::Ion, v::Vector)
    return map(x -> clear_sublevel_alias!(I, x), v)
end

"""
    clear_all_sublevel_aliases!(I::Ion)
Applies `clear_sublevel_alias!` to all sublevels of `I`.
"""
function clear_all_sublevel_aliases!(I::Ion)
    return empty!(I.sublevel_aliases)
end

"""
    set_stark_shift!(I::Ion, sublevel, shift::Real)
Applies a stark shift `shift` to the chosen `sublevel` of `I` (overwriting any previously assigned stark shift).
"""
function set_stark_shift!(I::Ion, sublevel::Tuple{String, Real}, shift::Real)
    validatesublevel(I, sublevel)
    return I.stark_shift[(sublevel[1], Rational(sublevel[2]))] = shift
end
set_stark_shift!(I::Ion, alias::String, shift::Real) =
    set_stark_shift!(I, alias2sublevel(I, alias), shift)
"""
    set_stark_shift!(I::Ion, stark_shift_dict::Dict)
Applies `set_stark_shift(I, sublevel, shift)` to all pairs `sublevel => shift` of the Dict `stark_shift_dict`.
"""
function set_stark_shift!(I::Ion, stark_shift_dict::Dict)
    for sublevel in keys(stark_shift_dict)
        set_stark_shift!(I, sublevel, stark_shift_dict[sublevel])
    end
end

"""
    zero_stark_shift!(I::Ion)
Sets the stark shift of all sublevels of `I` to zero.
"""
function zero_stark_shift!(I::Ion)
    for sublevel in keys(stark_shift(I))
        I.stark_shift[sublevel] = 0.0
    end
end

#############################################################################################
# Properties of ion electronic levels and sublevels
#############################################################################################

"""
    levels(I::Ion)
Returns array of all energy levels of `I`.
"""
levels(I::Ion) = unique([sublevel[1] for sublevel in sublevels(I)])

"""
    sublevelalias(I::Ion, sublevel::Tuple{String,Real})
Returns the alias assined to `sublevel` of `I`. If no alias is assigned, returns `nothing`.
"""
function sublevelalias(I::Ion, sublevel::Tuple{String, Real})
    validatesublevel(I, sublevel)
    alias_dict = sublevel_aliases(I)
    aliases = [k for (k, v) in alias_dict if v == sublevel]
    if length(aliases) == 0
        return nothing
    elseif length(aliases) == 1
        return aliases[1]
    else
        @warn "multiple aliases point to the same level $sublevel"
    end
end

"""
    alias2sublevel(I::Ion, alias::String)
Returns the sublevel corresponding to the given alias `alias` of `I`. Inverse function of `sublevelalias`.
"""
function alias2sublevel(I::Ion, alias::String)
    all_aliases = I.sublevel_aliases
    @assert alias in keys(all_aliases) "Ion does not contain any sublevel with the alias $alias. Use sublevel_aliases(Ion) to see available aliases."
    return all_aliases[alias]
end

"""
    sublevel2level(I::Ion, sublevel)
Retuns the energy level of `I` corresponding to `sublevel`.
"""
function sublevel2level(I::Ion, sublevel::Tuple{String, Real})
    validatesublevel(I, sublevel)
    return sublevel[1]
end
function sublevel2level(I::Ion, alias::String)
    validatesublevel(I, alias)
    sublevel = alias2sublevel(I, alias)
    return sublevel[1]
end

# quantumnumbers is written to be able to accept either a level or sublevel in the second argument
# Since both levels and aliases are strings, multidispatach can't tell the difference, so the second method distinguishes these cases with an if statement.
"""
    quantumnumbers(I::Ion, level::String)
    quantumnumbers(I::Ion, sublevel)
Returns the quantum numbers of an energy level or sublevel of `I` as a `NamedTuple`.

If second argument is a level, returns `(:n, :i, :s, :l, :j, :f)`

If second argument is a sublevel, returns `(:n, :i, :s, :l, :j, :f, :m)`
"""
function quantumnumbers(I::Ion, sublevel::Tuple{String, Real})
    validatesublevel(I, sublevel)
    levelstruct = speciesproperties(I).full_level_structure[sublevel[1]]
    names = (:n, :i, :s, :l, :j, :f, :m)
    values = (
        levelstruct.n,
        nuclearspin(I),
        1 // 2,
        levelstruct.l,
        levelstruct.j,
        levelstruct.f,
        Rational(sublevel[2])
    )
    return NamedTuple{names}(values)
end
function quantumnumbers(I::Ion, level_or_alias::String)
    # If the second argument is a String, it could be either a level name or the alias of a sublevel
    if level_or_alias in levels(I)
        # Second argument is a level name. Leave out the m quantum number
        levelstruct = speciesproperties(I).full_level_structure[level_or_alias]
        names = (:n, :i, :s, :l, :j, :f)
        values = (
            levelstruct.n,
            nuclearspin(I),
            1 // 2,
            levelstruct.l,
            levelstruct.j,
            levelstruct.f
        )
        return NamedTuple{names}(values)
    else
        # Second argument is a sublevel alias.
        quantumnumbers(I, alias2sublevel(I, level_or_alias))
    end
end

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
    stark_shift(I::Ion, sublevel)
Returns the assigned stark shift of `sublevel` of Ion `I`.
"""
function stark_shift(I::Ion, sublevel::Tuple{String, Real})
    validatesublevel(I, sublevel)
    return stark_shift(I)[sublevel]
end
stark_shift(I::Ion, alias::String) = stark_shift(I, alias2sublevel(I, alias))

"""
    zeeman_shift(I::Ion, sublevel, B::Real)
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
    energy(I::Ion, level::String)
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
    leveltransitions(I::Ion)
Returns all allowed transitions between levels of `I` as a vector of `Tuple{String,String}`.
"""
function leveltransitions(I::Ion)
    list = []
    lvls = levels(I)
    for levelpair in keys(speciesproperties(I).full_transitions)
        if levelpair[1] in lvls && levelpair[2] in lvls
            push!(list, levelpair)
        end
    end
    return list
end

"""
    subleveltransitions(I::Ion)
Returns all allowed transitions between sublevels of `I` as a vector of `Tuple{S,S}` where `S=Tuple{String,Real}`.
"""
function subleveltransitions(I::Ion)
    list = []
    for transition in leveltransitions(I)
        (L1, L2) = transition
        multipole = transitionmultipole(I, L1, L2)
        sublevels1 = [sublevel for sublevel in sublevels(I) if sublevel[1] == L1]
        sublevels2 = [sublevel for sublevel in sublevels(I) if sublevel[1] == L2]
        for sl1 in sublevels1
            for sl2 in sublevels2
                m1 = sl1[2]
                m2 = sl2[2]
                if abs(m2 - m1) <= parse(Int, multipole[2])
                    # Only add to list of sublevel transitions if Δm is not larger than the transition multipole allows (1 for E1, 2 for E2, etc)
                    push!(list, (sl1, sl2))
                end
            end
        end
    end
    return list
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

#############################################################################################
# Functions for constructing ion structs
#############################################################################################

function _construct_sublevels(selected_sublevels, properties)
    full_level_structure = properties.full_level_structure

    # If selected_sublevels is blank, use the default selection. If it is "all", use all sublevels.
    if selected_sublevels === nothing
        if !ismissing(properties.default_sublevel_selection)
            selected_sublevels = properties.default_sublevel_selection
        else
            @error "no level structure specified in constructor, and no default level structure specified for this ion species"
        end
    elseif selected_sublevels == "all"
        selected_sublevels = [(sublevel, "all") for sublevel in keys(full_level_structure)]
    end

    # Construct the list of sublevels
    sublevels = []
    for manifold in selected_sublevels
        # Ensure that the string is a valid level
        level = manifold[1]
        @assert level in keys(full_level_structure) "invalid level $level"
        @assert level ∉ [k[1] for k in keys(sublevels)] "multiple instances of level $level in ion constructor call"
        level_structure = full_level_structure[level]

        # Add chosen sublevels
        selectedms = manifold[2]
        f = level_structure.f

        m_allowed = Array((-f):f)
        if selectedms == "all"
            selectedms = m_allowed
        elseif !(typeof(selectedms) <: Array)
            selectedms = [selectedms]
        end
        for m in selectedms
            m = Rational(m)
            @assert m in m_allowed "Zeeman sublevel m = $m not allowed for state $level with f = $f"
            @assert (level, m) ∉ sublevels "repeated instance of sublevel $m in state $level"
            push!(sublevels, (level, m))
        end
    end

    return sublevels
end

function _construct_starkshift(starkshift, sublevels)
    starkshift_full = OrderedDict{Tuple, Real}()
    for sublevel in sublevels
        starkshift_full[sublevel] =
            (haskey(starkshift, sublevel) ? starkshift[sublevel] : 0.0)
    end
    return starkshift_full
end

#############################################################################################
# Overrides of Base functions
#############################################################################################

Base.getindex(I::Ion, state::Union{Tuple{String, Real}, String, Int}) = ionstate(I, state)

function Base.getproperty(I::Ion, s::Symbol)
    if s == :ionnumber || s == :position
        if typeof(getfield(I, s)) <: Missing
            @warn "ion has not been added to a configuration"
            return missing
        end
    end
    return getfield(I, s)
end

function Base.setproperty!(I::Ion, s::Symbol, v::Tv) where {Tv}
    if (s == :species_properties || s == :shape || s == :number || s == :position)
        return
    elseif s == :sublevels
        Core.setproperty!(I, :sublevels, _construct_sublevels(v, speciesproperties(I)))
        # Also update the stark shift dict as necessary; keep old stark shift values and assign zero stark shift to new sublevels
        starkshift_full_old = stark_shift(I)
        starkshift_full_new = OrderedDict{Tuple, Real}()
        for sublevel in sublevels(I)
            starkshift_full_new[sublevel] = (
                haskey(starkshift_full_old, sublevel) ? starkshift_full_old[sublevel] : 0.0
            )
        end
        Core.setproperty!(I, :stark_shift, starkshift_full_new)
        return
    elseif s == :sublevel_aliases
        for (alias, sublevel) in v
            validatesublevel(I, sublevel)
        end
        Core.setproperty!(I, :sublevel_aliases, v)
    elseif s == :stark_shift
        starkshift_full_new = stark_shift(I) # New stark shift dict is initially identical to the old one
        for (sublevel, shift) in v
            validatesublevel(I, sublevel)
            starkshift_full_new[sublevel] = shift # Change the shifts as necessary
        end
        Core.setproperty!(I, :stark_shift, starkshift_full_new)
    end
end

function Base.:(==)(b1::T, b2::T) where {T <: Ion}
    # Takes two ions to be equal if they are the same species, contain the same sublevels, and have the same stark shifts
    # Does not care about sublevel aliases or the ordering of sublevels
    return (
        b1.species_properties == b2.species_properties &&
        sort(b1.sublevels) == sort(b2.sublevels) &&
        sort(b1.stark_shift) == sort(b2.stark_shift)
    )
end

"""
IonProperties type.
**Required keywords**
* `mass`: Mass of ion in kg
* `charge`: Charge of ion in units of elementary charge
* `full_level_structure`: OrderedDict describing properties of each energy level
  * `key::String`: Name of energy level. Spectroscopic notation encouraged, e.g. `"S1/2,f=1"`
  * `value::NamedTuple(:n, :l, :j, :f, :E)`: Quantum numbers `n`, `l`, `j`, `f`, and energy `E` (in Hz)
* `full_transitions`: Dict of all allowed transitions between energy levels
  * `key::Tuple{String,String}` Pair of levels, ordered (lower, upper) in energy
  * `value::NamedTuple(:multipole, :einsteinA)`: Leading-order multipole of the transition (e.g. `"E1"`, `"E2"`) and Einstein A coefficient (between fine structure levels only; hyperfine factors are calculated when needed)

 **Optional keywords**
 * `default_sublevel_selection`: Default value of `selected_sublevels` argument in Ion constructor
 * `gfactors`: `Dict(level::String => g::Real)` Custom Landé g-factors, if contributions from higher-than-first-order perturbations are desired
 * `nonlinear_zeeman`: `Dict` describing nonlinear contributions to Zeeman shift of certain sublevels
   * `key::Tuple{String,Real}`: sublevel name
   * `value::Function(B::Real)`: Nonlinear term(s) of Zeeman shift. Full Zeeman shift will be calculated as the sum of the usual linear term and this function
"""
struct IonProperties
    shortname::String
    mass::Number
    charge::Number
    nuclearspin::Number
    full_level_structure::OrderedDict{
        String,
        NamedTuple{(:n, :l, :j, :f, :E), T} where T <: Tuple
    }
    full_transitions::Dict{
        Tuple{String, String},
        NamedTuple{(:multipole, :einsteinA), Tuple{String, Float64}}
    }
    default_sublevel_selection::Union{Vector{Tuple{String, String}}, Missing}
    gfactors::Union{Dict{String, Number}, Missing}
    nonlinear_zeeman::Union{Dict{Tuple{String, Real}, Function}, Missing}
end
IonProperties(;
    shortname,
    mass,
    charge,
    nuclearspin,
    full_level_structure,
    full_transitions,
    default_sublevel_selection = missing,
    gfactors = missing,
    nonlinear_zeeman = missing
) = IonProperties(
    shortname,
    mass,
    charge,
    nuclearspin,
    full_level_structure,
    full_transitions,
    default_sublevel_selection,
    gfactors,
    nonlinear_zeeman
)

function Base.:(==)(p1::IonProperties, p2::IonProperties)
    for fn in fieldnames(IonProperties)
        p1fn = getfield(p1, fn)
        p2fn = getfield(p2, fn)
        if ismissing(p1fn) && ismissing(p2fn)
            continue
        elseif p1fn == p2fn
            continue
        else
            return false
        end
    end
    return true
end

"""
IonInstance(selected_sublevels::Vector{Tuple}[, starkshift::Dict])
Ion instance of some species

`selected_sublevels` specifies which energy sublevels will be present in the Hilbert space of this Ion instance, as a subset of all possible sublevels.

Each element of `selected_sublevels` is a 2-element Tuple (level, sublevels), with the first element being the name of a level and the second specifying which sublevels should be included.
Allowed sublevels are those whose magnetic quantum number `m` is in the set {`-f`, `-f+1`, `-f+2`, ... `f-1`, `f`}, where `f` is the total angular momentum quantum number of `level`.
For each `level` specified there are three allowed options to specify the set of `sublevels` to include:
* `sublevels::Real`: Includes only one `m = sublevels`
* `sublevels::Vector{Real}`: Includes all sublevels whose magnetic quantum number `m` is in `sublevels`
* `sublevels = "all"`: Includes all allowed sublevels

If instead `selected_sublevels = "all"`, then all sublevels of all levels are included.

Omission of a level in `selected_sublevels` will exclude all sublevels.

**Fields**
* `species_properties::NamedTuple`: Contains constants specifying parameters specific to species
* `sublevels`::Vector{Tuple{String,Real}}: List of all sublevels present in the Hilbert space
* `sublevel_aliases::Dict{String,Tuple}`: Dict specifying aliases assigned to sublevels, in the format `alias => sublevel`
* `shape`::Vector{Int}: Dimension of the Hilbert space
* `stark_shift::OrderedDict`: A dictionary with keys denoting the selected levels and values, a real number for describing a shift of the level's energy. This is just a convenient way to add Stark shifts to the simulation without additional resources
* `ionnumber`: When the ion is added to an `IonConfiguration`, this value keeps track of its order in the chain
* `position`: When the ion is added to an `IonConfiguration`, this value keeps track of its physical position in meters
"""
mutable struct IonInstance{Species <: Any} <: Ion
    # fields
    species_properties::IonProperties
    sublevels::Vector{Tuple{String, Real}}
    sublevel_aliases::Dict{String, Tuple}
    shape::Vector{Int}
    stark_shift::OrderedDict{Tuple, Real}
    ionnumber::Union{Int, Missing}
    position::Union{Real, Missing}

    # constructors (overrides default)
    function IonInstance{Species}(
        properties,
        selected_sublevels::Union{Vector{Tuple{String, T}}, String, Nothing} where {T} = nothing,
        starkshift = Dict()
    ) where {Species <: Any}
        sublevels = _construct_sublevels(selected_sublevels, properties)
        shape = [length(sublevels)]
        starkshift_full = _construct_starkshift(starkshift, sublevels)
        return new{Species}(
            properties,
            sublevels,
            Dict(),
            shape,
            starkshift_full,
            missing,
            missing
        )
    end
    function IonInstance{Species}(
        species_properties,
        sublevels,
        sublevel_aliases,
        shape,
        stark_shift,
        ionnumber,
        position
    ) where {Species <: Any}
        sublevels = deepcopy(sublevels)
        sublevel_aliases = deepcopy(sublevel_aliases)
        shape = copy(shape)
        stark_shift = deepcopy(stark_shift)
        return new{Species}(
            species_properties,
            sublevels,
            sublevel_aliases,
            shape,
            stark_shift,
            ionnumber,
            position
        )
    end
end

function Base.print(I::IonInstance)
    println(I.species_properties.shortname + "\n")
    for (k, v) in I.selected_sublevel_structure
        println(k, ": ", v)
    end
end

Base.show(io::IO, I::IonInstance) = println(io, I.species_properties.shortname)  # suppress long output
