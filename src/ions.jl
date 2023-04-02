using WignerSymbols: wigner3j, wigner6j
using LinearAlgebra: cross
using IonSim.PhysicalConstants

export Ion,
    # SpeciesProperties,
    speciesproperties,
    sublevels,
    sublevelaliases,
    sublevelalias,
    shape,
    manualshift,
    ionnumber,
    ionposition,
    mass,
    charge,
    nuclearspin,
    zeromanualshift!,
    manualshift!,
    sublevel,
    level,
    sublevelalias!,
    clearsublevelalias!,
    levels,
    quantumnumbers,
    landegf,
    zeemanshift,
    energy,
    transitionfrequency,
    transitionwavelength,
    leveltransitions,
    subleveltransitions,
    einsteinA,
    transitionmultipole,
    lifetime,
    matrixelement,
    IonInstance


#############################################################################################
# Ion - the physical parameters defining an ion's structure
#############################################################################################

"""
    Ion

Contains the physical parameters defining the internal structure of an isolated ion.
"""
abstract type Ion <: IonSimBasis end


#############################################################################################
# Ion fields
#############################################################################################

"""
    speciesproperties(ion::Ion)

Returns the `ion.speciesproperties`, which contains species-specific information including 
mass, energy levels, and transitions.
"""
speciesproperties(I::Ion) = I.speciesproperties

"""
    sublevels(ion::Ion)

Returns the energy sublevels in the Hilbert space of `ion`. Equivalent to `ion.sublevels`.
"""
sublevels(I::Ion) = I.sublevels

"""
    sublevelaliases(ion::Ion)

Returns a `Dict` specifying all aliases assigned to sublevels of `ion`, in the format
`alias => sublevel`, `ion.sublevelalias`
"""
sublevelaliases(I::Ion) = I.sublevelaliases

"""
    shape(ion::Ion)

Returns the dimension of the ion's Hilbert space `ion.shape`
"""
shape(I::Ion) = I.shape

"""
    manualshift(ion::Ion)

Returns the Dict of manualshift for `ion`'s energy levels. Equivalent to `ion.manualshift`.
"""
manualshift(I::Ion) = I.manualshift

# Will be removing this, instead the ion will have a field pointing to a particular IonTrap
"""
    ionnumber(ion::Ion)

Returns `ion`'s number in its IonTrap `ion.ionnumber`
If `ion` has not been added to an IonTrap, returns `missing`.
"""
function ionnumber(I::Ion)
    if typeof(I.ionnumber) <: Missing
        @warn "ion has not been added to an IonTrap"
        return missing
    else
        return I.ionnumber
    end
end

# Will also be removing this 
"""
    ionposition(ion::Ion)

Returns `ion`'s position in its IonTrap, `ion.ionposition` in m.
If `ion` has not been added to an IonTrap, returns `missing`.
"""
function ionposition(I::Ion)
    if typeof(I.ionposition) <: Missing
        @warn "ion has not been added to an IonTrap"
        return missing
    else
        return I.ionposition
    end
end

"""
    mass(ion::Ion)

Returns the mass of `ion` in kg. Equivalent to `ion.speciesproperties.mass`.
"""
mass(I::Ion) = speciesproperties(I).mass

"""
    charge(ion::Ion)

Returns the charge of `ion` in Coulombs. 
Equivalent to `ion.speciesproperties.charge * PhysicalConstants.e`.
"""
charge(I::Ion) = speciesproperties(I).charge * e

"""
    nuclearspin(ion::Ion)

Returns the nuclear spin of `ion`. Equivalent to `ion.speciesproperties.nuclearspin`.
"""
nuclearspin(I::Ion) = speciesproperties(I).nuclearspin


#############################################################################################
# SpeciesProperties - A container for the intrinsic physical properties of a particular 
#                     ion species.
#############################################################################################

struct SpeciesProperties
    shortname::String  # an intuitive string identifier for the ion
    mass::Number  # mass of ion in kg
    charge::Number  # charge of ion in units of elementary charge e
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

function SpeciesProperties(;
        shortname,
        mass,
        charge,
        nuclearspin,
        full_level_structure,
        full_transitions,
        default_sublevel_selection=missing,
        gfactors=missing,
        nonlinear_zeeman=missing
    )  
    SpeciesProperties(
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
end


#############################################################################################
# IonInstance
#############################################################################################

"""
    IonInstance(selected_sublevels::Vector{Tuple}[, manualshift::Dict])

Ion instance of some species.

`selected_sublevels` specifies which energy sublevels will be present in the Hilbert space of 
this Ion instance, as a subset of all possible sublevels.

Each element of `selected_sublevels` is a 2-element Tuple `(level::String, sublevels)`, with 
the first element being the name of a level and the second specifying which sublevels should 
be included. Allowed sublevels are those whose magnetic quantum number `m` is in the set 
{`-f`, `-f+1`, `-f+2`, ... `f-1`, `f`}, where `f` is the total angular momentum quantum 
number of `level`. For each `level` specified there are three allowed options to specify the 
set of `sublevels` to include:

* `sublevels::Real`: Includes only one `m = sublevels`
* `sublevels::Vector{Real}`: Includes all sublevels whose magnetic quantum number `m` is in 
    `sublevels`
* `sublevels = "all"`: Includes all allowed sublevels
If an element of `selected_sublevels` utilizes the first option, specifying a single `m`, one 
may optionally may this a 3-element tuple instead: `(level::String, m::Real, alias::String)`, 
assinging this particular sublevel the alias `alias`.

Omission of a level in `selected_sublevels` will exclude all sublevels.

**Fields**
* `speciesproperties`: A container of species-specific properties
* `sublevels`::Vector{Tuple{String,Real}}: List of all sublevels present in the Hilbert space
* `sublevelaliases::Dict{String,Tuple}`: Dict specifying aliases assigned to sublevels, in the 
    format `alias => sublevel`
* `shape`::Vector{Int}: Dimension of the Hilbert space
* `manualshift::OrderedDict`: A dictionary with keys denoting the selected levels and values, 
    a real number for describing a shift of the level's energy. This is just a convenient way 
    to add manual shifts to the simulation, such as Stark shifts off of energy levels 
    not present in the Hilbert space, without additional resources
* `ionnumber`: When the ion is added to an `IonTrap`, this value keeps track of its order
* `ionposition`: When the ion is added to an `IonTrap`, this value keeps track of its physical 
    position in meters
"""
mutable struct IonInstance{Species<:Any} <: Ion
    # fields
    speciesproperties::SpeciesProperties
    sublevels::Vector{Tuple{String, Real}}
    sublevelaliases::Dict{String, Tuple}
    shape::Vector{Int}
    manualshift::OrderedDict{Tuple, Real}
    ionnumber::Union{Int, Missing}
    ionposition::Union{Real, Missing}

    # constructors (overrides default)
    function IonInstance{Species}(
        properties,
        selected_sublevels=nothing,
        manualshift=Dict()
    ) where {Species<:Any}
        (sublevels, aliases) = _construct_sublevels(selected_sublevels, properties)
        shape = [length(sublevels)]
        manualshift_full = _construct_manualshift(manualshift, sublevels)
        return new{Species}(
            properties,
            sublevels,
            aliases,
            shape,
            manualshift_full,
            missing,
            missing
        )
    end

    function IonInstance{Species}(
        speciesproperties,
        sublevels,
        sublevelaliases,
        shape,
        manualshift,
        ionnumber,
        ionposition
    ) where {Species<:Any}
        sublevels = deepcopy(sublevels)
        sublevelaliases = deepcopy(sublevelaliases)
        shape = copy(shape)
        manualshift = deepcopy(manualshift)
        return new{Species}(
            speciesproperties,
            sublevels,
            sublevelaliases,
            shape,
            manualshift,
            ionnumber,
            ionposition
        )
    end
end

#############################################################################################
# Functions for constructing ion structs
#############################################################################################

function _construct_sublevels(selected_sublevels, properties)
    full_level_structure = properties.full_level_structure

    # If selected_sublevels is blank, use the default selection. If it is "all", use all 
    # sublevels.
    if isnothing(selected_sublevels)
        if !ismissing(properties.default_sublevel_selection)
            selected_sublevels = properties.default_sublevel_selection
        else
            @error (
                "no level structure specified in constructor, and no default level structure \
                specified for this ion species"
            )
        end
    elseif selected_sublevels == "all"
        selected_sublevels = [(sublevel, "all") for sublevel in keys(full_level_structure)]
    end

    # Construct the list of sublevels
    sublevels = []
    aliases = Dict()
    for manifold in selected_sublevels
        # Ensure that the string is a valid level
        level = manifold[1]
        @assert level in keys(full_level_structure) "invalid level $level"
        level_structure = full_level_structure[level]

        # Create aliases, if applicable
        selectedms = manifold[2]
        if length(manifold) == 3
            @assert typeof(selectedms) <: Real (
                    "$manifold: sublevel aliases may only be assigned to a single sublevel"
            )
            sl = (level, Rational(selectedms))
            alias = manifold[3]
            aliases[alias] = sl
        end

        # Add chosen sublevels
        f = level_structure.f

        m_allowed = Array((-f):f)
        if selectedms == "all"
            selectedms = m_allowed
        elseif !(typeof(selectedms) <: Array)
            selectedms = [selectedms]
        end
        for m in selectedms
            m = Rational(m)
            @assert m in m_allowed (
                "Zeeman sublevel m = $m not allowed for state $level with f = $f"
            )
            @assert (level, m) ∉ sublevels "repeated instance of sublevel $m in state $level"
            push!(sublevels, (level, m))
        end
    end

    return (sublevels, aliases)
end

function _construct_manualshift(manualshift, sublevels)
    manualshift_full = OrderedDict{Tuple, Real}()
    for sublevel in sublevels
        manualshift_full[sublevel] =
            (haskey(manualshift, sublevel) ? manualshift[sublevel] : 0.0)
    end
    return manualshift_full
end


#############################################################################################
# Functions to modify ion fields
#############################################################################################

validatesublevel(I::Ion, sublevel::Tuple{String, Real}) =
    @assert sublevel in sublevels(I) (
            "Ion does not contain sublevel $sublevel. Use sublevels(Ion) to see list of 
            available sublevels."
        )
validatesublevel(I::Ion, alias::String) = validatesublevel(I, sublevel(I, alias))

"""
    sublevelalias!(I::Ion, sublevel::Tuple{String,Real}, alias::String)

Assigns an alias `alias` to `sublevel` of `I`. This then allows one to pass `alias` in place 
of `sublevel` (for `I` only) into any function which accepts a sublevel as an argument.
"""
function sublevelalias!(I::Ion, sublevel::Tuple{String, Real}, alias::String)
    validatesublevel(I, sublevel)
    @assert alias ∉ levels(I) "Cannot make alias name identical to level name ($alias)"
    sublevel_rational = (sublevel[1], Rational(sublevel[2]))   # Force m to be Rational
    return I.sublevelaliases[alias] = sublevel_rational
end

"""
    sublevelalias!(I::Ion, aliasassignments)

`aliasassignments` specifies which aliases should be assigned to which sublevels. There are 
two options to do this:
* `aliasassignments` is a Vector of Tuples, with the first element of each being the sublevel 
 (`sublevel::Tuple{String,Real}`) and the second being its assigned alias (`alias::String`)
* `aliasassignments` is a Dict with the format `alias::String => sublevel::Tuple{String,Real}`
Calls `sublevelalias!(I, sublevel, alias)` for each pair `sublevel, alias`.
"""
function sublevelalias!(
    I::Ion,
    pairs::Vector{Tuple{Tuple{String, R}, String}} where {R <: Real}
)
    for (sublevel, alias) in pairs
        sublevelalias!(I, sublevel, alias)
    end
end

function sublevelalias!(I::Ion, aliasdict::Dict{String, Tuple{String, R}} where {R <: Real})
    for (alias, sublevel) in aliasdict
        sublevelalias!(I, sublevel, alias)
    end
end

"""
    clearsublevelalias!(I::Ion, [, sublevel])

Erases the assignment of an alias to `sublevel` of Ion `I`. Accepts either the full sublevel 
`Tuple{String,Real}` or its alias `String`.
Also accepts a vector of sublevels to clear multiple alias assignments in a single call.
"""
function clearsublevelalias!(I::Ion, sublevel::Tuple{String, Real})
    alias = sublevelalias(I, sublevel)
    return delete!(I.sublevelaliases, alias)
end
clearsublevelalias!(I::Ion, alias::String) = delete!(I.sublevelaliases, alias)
clearsublevelalias!(I::Ion, v::Vector) = map(x -> clearsublevelalias!(I, x), v)

"""
If sublevel not specified, the assignment of all sublevel aliases of Ion `I` are erased. 
"""
clearsublevelalias!(I::Ion) = return empty!(I.sublevelaliases)

"""
    manualshift!(I::Ion, sublevel, shift::Real)

Applies a manual shift `shift` to the chosen `sublevel` of `I` (overwriting any previously 
assigned manual shift).
"""
function manualshift!(I::Ion, sublevel::Tuple{String, Real}, shift::Real)
    validatesublevel(I, sublevel)
    return I.manualshift[(sublevel[1], Rational(sublevel[2]))] = shift
end
manualshift!(I::Ion, alias::String, shift::Real) = manualshift!(I, sublevel(I, alias), shift)

"""
    manualshift!(I::Ion, manualshift_dict::Dict)

Applies `manualshift(I, sublevel, shift)` to all pairs `sublevel => shift` of the Dict 
`manual_shift_dict`.
"""
function manualshift!(I::Ion, manual_shift_dict::Dict)
    for sublevel in keys(manual_shift_dict)
        manualshift!(I, sublevel, manual_shift_dict[sublevel])
    end
end

"""
    zeromanualshift!(I::Ion)

Sets the manual shift of all sublevels of `I` to zero.
"""
function zeromanualshift!(I::Ion)
    for sublevel in keys(manualshift(I))
        manualshift!(I, sublevel, 0.0)
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
    alias_dict = sublevelaliases(I)
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
    sublevel(I::Ion, alias::String)

Returns the sublevel corresponding to the given alias `alias` of `I`. 
Inverse function of `sublevelalias`.
"""
function sublevel(I::Ion, alias::String)
    all_aliases = I.sublevelaliases
    @assert alias in keys(all_aliases) (
            "Ion does not contain any sublevel with the alias $alias. Use sublevelaliases(Ion) 
            to see available aliases."
    )
    return all_aliases[alias]
end

"""
    level(I::Ion, sublevel)

Retuns the energy level of `I` corresponding to `sublevel`.
"""
function level(I::Ion, sublevel::Tuple{String, Real})
    validatesublevel(I, sublevel)
    return sublevel[1]
end

function level(I::Ion, alias::String)
    validatesublevel(I, alias)
    sl = sublevel(I, alias)
    return sl[1]
end

#= 
quantumnumbers is written to be able to accept either a level or sublevel in the second 
argument. Since both levels and aliases are strings, multidispatach can't tell the difference, 
so the second method distinguishes these cases with an if statement.
=#
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
    # If the second argument is a String, it's either a level name or the alias of a sublevel
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
        quantumnumbers(I, sublevel(I, level_or_alias))
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
landegj(l::Real, j::Real, s::Real=1 // 2) =
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
landegf(l::Real, j::Real, f::Real, i::Real, s::Real=1 // 2) =
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
    manualshift(I::Ion, sublevel)

Returns the assigned manual shift of `sublevel` of Ion `I`.
"""
function manualshift(I::Ion, sublevel::Tuple{String, Real})
    validatesublevel(I, sublevel)
    return manualshift(I)[sublevel]
end
manualshift(I::Ion, alias::String) = manualshift(I, sublevel(I, alias))

"""
    zeemanshift(I::Ion, sublevel, B::Real)

Returns the Zeeman shift at a magnetic field of `B` of `sublevel` of `I`.
If `sublevel` has a custom g-factor defined, then this is used. Otherwise, `landegf` is used 
to compute the Landé g-factor. Zeeman shift calculated as ``ΔE = (μ_B/ħ) ⋅ g_f ⋅ B ⋅ m / 2π``
"""
function zeemanshift(I::Ion, sublevel::Tuple{String, Real}, B::Real)
    validatesublevel(I, sublevel)
    properties = speciesproperties(I)
    if !ismissing(properties.nonlinear_zeeman) &&
       haskey(properties.nonlinear_zeeman, sublevel)
        nonlinear = properties.nonlinear_zeeman[sublevel](B)
    else
        nonlinear = 0.0
    end
    return zeemanshift(B, landegf(I, sublevel[1]), sublevel[2]) + nonlinear
end
zeemanshift(B::Real, g::Real, m::Real) = (μB / ħ) * g * B * m / 2π

zeemanshift(B::Real, l::Real, j::Real, f::Real, m::Real, i::Real, s::Real=1 // 2) =
    zeemanshift(B, landegf(l, j, f, i, s), m)

zeemanshift(B::Real, qnums::NamedTuple) =
    zeemanshift(B, qnums.l, qnums.j, qnums.f, qnums.m, qnums.i, qnums.s)

zeemanshift(I::Ion, alias::String, B::Real) = zeemanshift(I, sublevel(I, alias), B)

#= 
This function is written to be able to accept either a level or sublevel in the second 
argument. Since both levels and aliases are strings, multidispatach can't tell the difference, 
so the second method distinguishes these cases with an if statement.
=#
"""
    energy(I::Ion, sublevel; B=0, ignore_manualshift=false)

Calculates the energy of `sublevel` of `I`. A Zeeman shift can be included by setting a nonzero 
value of `B`. The manual shift can be omitted by setting `ignore_manualshift=true`.
"""
function energy(I::Ion, sublevel::Tuple{String, Real}; B=0, ignore_manualshift=false)
    validatesublevel(I, sublevel)
    E0 = speciesproperties(I).full_level_structure[sublevel[1]].E
    zeeman = zeemanshift(I, sublevel, B)
    manual = (ignore_manualshift ? 0.0 : manualshift(I, sublevel))
    return E0 + zeeman + manual
end

"""
    energy(I::Ion, level::String)

Returns the energy of `level` of `I`.
"""
function energy(I::Ion, level_or_alias::String; B=0, ignore_manualshift=false)
    # If the second argument is a String, it could be either a level name or the alias of a 
    # sublevel
    if level_or_alias in levels(I)
        # Second argument is a level name. Return the bare energy of that level.
        return speciesproperties(I).full_level_structure[level_or_alias].E
    else
        # Second argument is a sublevel alias.
        return energy(
            I,
            sublevel(I, level_or_alias),
            B=B,
            ignore_manualshift=ignore_manualshift
        )
    end
end

"""
    transitionfrequency(I::Ion, transition::Tuple; B=0, ignore_manualshift=false)

`transition` is a Tuple of two sublevels or levels. Computes the absolute values of the 
difference in energies between `transition[1]` and `transition[2]`. If between sublevels, 
then the Zeeman shift may be included by setting the value of the magnetic field `B`, and 
manual shifts may be omitted by setting `ignore_manualshift=true`.
"""
function transitionfrequency(I::Ion, transition::Tuple; B=0, ignore_manualshift=false)
    # Multidispatch of the function energy should make this work regardless of whether the 
    # transition is between levels or sublevels, and regardless of whether or not aliases are 
    # used
    return abs(
        energy(I, transition[1], B=B, ignore_manualshift=ignore_manualshift) -
        energy(I, transition[2], B=B, ignore_manualshift=ignore_manualshift)
    )
end

"""
    transitionwavelength(I::Ion, transition::Tuple; B=0, ignore_manualshift=false)

Returns the wavelength corresponding to `transitionfrequency(I::Ion, transition::Tuple; B=0, 
ignore_manualshift=false)`.
"""
function transitionwavelength(I::Ion, transition::Tuple; B=0, ignore_manualshift=false)
    return c /
           transitionfrequency(I, transition, B=B, ignore_manualshift=ignore_manualshift)
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

Returns all allowed transitions between sublevels of `I` as a vector of `Tuple{S,S}` where 
`S=Tuple{String,Real}`.
"""
function subleveltransitions(I::Ion)
    list = []
    for transition in leveltransitions(I)
        (L1, L2) = transition
        multipole = transitionmultipole(I, (L1, L2))
        sublevels1 = [sublevel for sublevel in sublevels(I) if sublevel[1] == L1]
        sublevels2 = [sublevel for sublevel in sublevels(I) if sublevel[1] == L2]
        for sl1 in sublevels1
            for sl2 in sublevels2
                m1 = sl1[2]
                m2 = sl2[2]
                if abs(m2 - m1) <= parse(Int, multipole[2])
                    #= Only add to list of sublevel transitions if Δm is not larger than the 
                       transition multipole allows (1 for E1, 2 for E2, etc)
                    =#
                    push!(list, (sl1, sl2))
                end
            end
        end
    end
    return list
end

"""
    einsteinA(I::Ion, leveltransition::Tuple)

Returns Einstein A coefficient corresponding to the transition 
`leveltransition[1] -> leveltransition[2]`. The first level must be the lower level and the 
second must be the upper level.
"""
function einsteinA(I::Ion, leveltransition::Tuple{String, String})
    (L1, L2) = leveltransition
    @assert (L1, L2) in leveltransitions(I) "invalid transition $L1 -> $L2"
    return speciesproperties(I).full_transitions[(L1, L2)].einsteinA
end

"""
    transitionmultipole(I::Ion, leveltransition::Tuple)

Returns the transition multiple (`'E1'`, `'E2'`, etc.) corresponding to the transition 
`leveltransition[1] -> leveltransition[2]`. The first level must be the lower level and the 
second must be the upper level.
"""
function transitionmultipole(I::Ion, leveltransition::Tuple{String, String})
    (L1, L2) = leveltransition
    @assert (L1, L2) in leveltransitions(I) "invalid transition $L1 -> $L2"
    return speciesproperties(I).full_transitions[(L1, L2)].multipole
end

"""
    lifetime(I::Ion, level::String)

Computes lifetime of `level` by summing the transition rates out of `level`. The sum is taken 
 over all levels that the species may have, rather than the levels present in the instance `I`.
"""
function lifetime(I::Ion, level::String)
    @assert level in keys(speciesproperties(I).full_level_structure) (
        "Ion species $(typeof(I)) does not contain level $level"
    )
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
    matrixelement(
        ion::Ion, transition::Tuple, I::Real, ϵhat::NamedTuple, 
        khat::NamedTuple, Bhat::NamedTuple=(;z=1)
    )
Computes the matrix elements (units of Hz) between two energy sublevels.

**args**
* `ion`: Ion undergoing transition
* `transition`: Tuple of sublevels (full names or aliases) between which the transition is 
            being calculated. Must be formatted such that 
            `energy(transition[2]) > energy(transition[1])`
* `I`: Intensity of the driving field
* `ϵhat`: Unit vector of light polarization
* `khat`: Unit vector of light wavevector
* `Bhat`: Unit vector of magnetic field
"""
function matrixelement(
        j1::Real,
        j2::Real,
        f1::Real,
        f2::Real,
        m1::Real,
        m2::Real,
        i::Real,
        ΔE::Real,
        A12::Real,
        multipole::String,
        I::Real,
        ϵhat::NamedTuple,
        khat::NamedTuple,
        Bhat::NamedTuple=(; z=1)
    )
    # Level 1 *must* be the lower level and level 2 *must* be the upper level
    # Note that in this function, i is the nuclear spin

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
        # Rotation matrix in axis-angle representation (axis=a, angle=theta)
        R = eye3 + sin(theta) * amatrix + (1 - cos(theta)) * amatrix^2
    end
    ϵhat_rotated = R * ϵhat_array
    khat_rotated = R * khat_array

    if multipole == "E1"
        if abs(q) > 1
            return 0
        else
            E = efield(I)
            hyperfine_factor =
                abs(sqrt((2 * f1 + 1) * (2 * f2 + 1)) * wigner6j(j2, i, f2, f1, 1, j1))
            geometric_factor = abs(
                sqrt(2j2 + 1) *
                wigner3j(f2, 1, f1, -m2, q, m1) *
                (transpose(c_rank1[q+2, :]) * ϵhat_rotated)
            )
            units_factor = abs(e * E / (2ħ) * sqrt(3 * A12 / (α * c * k^3)))
            return units_factor * hyperfine_factor * geometric_factor / 2π
        end
    elseif multipole == "E2"
        if abs(q) > 2
            return 0
        else
            E = efield(I)
            hyperfine_factor =
                abs(sqrt((2 * f1 + 1) * (2 * f2 + 1)) * wigner6j(j2, i, f2, f1, 2, j1))
            geometric_factor = abs(
                sqrt(2j2 + 1) *
                wigner3j(f2, 2, f1, -m2, q, m1) *
                (transpose(khat_rotated) * c_rank2[:, :, q+3] * ϵhat_rotated)
            )
            units_factor = abs(e * E / (2ħ) * sqrt(15 * A12 / (α * c * k^3)))
            return units_factor * hyperfine_factor * geometric_factor / 2π
        end
    else
        @error (
            "calculation of atomic transition matrix element for transition type $multipole \
            not currently supported"
        )
    end
end

function matrixelement(
        ion::Ion,
        transition::Tuple,
        I::Real,
        ϵhat::NamedTuple,
        khat::NamedTuple,
        Bhat::NamedTuple=(; z=1)
    )
    SL1 = transition[1]
    SL2 = transition[2]
    L1 = level(ion, SL1)
    L2 = level(ion, SL2)

    E1 = energy(ion, L1)
    E2 = energy(ion, L2)
    @assert E2 > E1 "transition must be formatted (lower level, upper level)"
    qn1 = quantumnumbers(ion, SL1)
    qn2 = quantumnumbers(ion, SL2)
    A12 = einsteinA(ion, (L1, L2))
    multipole = transitionmultipole(ion, (L1, L2))

    return matrixelement(
        qn1.j,
        qn2.j,
        qn1.f,
        qn2.f,
        qn1.m,
        qn2.m,
        nuclearspin(ion),
        E2 - E1,
        A12,
        multipole,
        I,
        ϵhat,
        khat,
        Bhat
    )
end


#############################################################################################
# Overrides of Base functions
#############################################################################################

Base.getindex(I::Ion, state::Union{Tuple{String, Real}, String, Int}) = ionstate(I, state)

function Base.print(ion::IonInstance)
    println(ion.speciesproperties.shortname)
    for sublevel in sublevels(ion)
        if sublevel in values(sublevelaliases(ion))
            alias = sublevelalias(ion, sublevel)
            println("$sublevel: \"$alias\"")
        else
            println(sublevel)
        end
    end
end

Base.show(io::IO, I::IonInstance) = print(io, I.speciesproperties.shortname)  # suppress long output

function Base.:(==)(p1::SpeciesProperties, p2::SpeciesProperties)
    for fn in fieldnames(SpeciesProperties)
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

function Base.:(==)(C1::Ion, C2::Ion)
   return (mass(C1) == mass(C2)) && (C1.speciesproperties == C2.speciesproperties)
end
