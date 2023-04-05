using WignerSymbols: wigner3j, wigner6j
using LinearAlgebra: cross
using IonSim.PhysicalConstants

export Ion,
    speciesproperties,
    selected_sublevels,
    shape,
    manualshift,
    ionnumber,  #
    ionposition, #
    mass,
    charge,
    nuclearspin,
    levels,
    availabletransitions, 
    IonInstance,
    einsteinA,
    transitionmultipole,
    subleveltransitions, 
    lifetime,
    manualshift!,
    zeromanualshift!,
    landegf,
    landegj,
    zeemanshift,
    energy,
    transitionfrequency, 
    transitionwavelength, 
    matrixelement


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
    selected_sublevels(ion::Ion)

Returns the selected energy sublevels in the Hilbert space of `ion`. Equivalent to 
`ion.sublevels`.
"""
selected_sublevels(I::Ion) = I.selected_sublevels

"""
    shape(ion::Ion)

Returns the dimension of the ion's Hilbert space. Equivalent to `ion.shape`
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
    nuclearspin::Number  # nuclear spin of the ion
    level_structure::OrderedDict{EnergyLevel, Real}  # Mapping of EnergyLevel to energy (scalar)
    transitions::Dict{  # Dictionary of allowed transitions between energy levels
        Tuple{EnergyLevel, EnergyLevel},
        NamedTuple{(:multipole, :einsteinA), Tuple{String, Float64}}
    }
    gfactors::Union{Dict{EnergyLevel, Number}, Missing}
    nonlinear_zeeman::Union{Dict{Tuple{EnergyLevel, Real}, Function}, Missing}
end

function SpeciesProperties(;
        shortname,
        mass,
        charge,
        nuclearspin,
        level_structure,
        transitions,
        gfactors=missing,
        nonlinear_zeeman=missing
    )
    hf_transitions = Dict()
    for (k, v) in transitions
        transition_info = v
        e1, e2 = k
        j1, j2 = e1.j, e2.j
        # l1, l2 = e1.l, e2.l

        multipole = v.multipole
        if multipole == "user-defined"
            @warn(
                "NotImplemented: User defined transition matrix elements are not yet \
                implemented. These will be skipped"
            )
            continue
        elseif !(multipole in ["E1", "E2", "M1"])
            @warn(
                "NotImplemented: Only E1, E2 and M1 transitions currently supported"
            )
            continue
        end
        if nuclearspin == 0
            hf_transitions = transitions
        else
            # F ∈ {|J-I|, |J-I| + 1, |J-I| + 2, ..., J+I}
            possible_F1 = [i for i in abs(j1 - nuclearspin):j1 + nuclearspin]
            possible_F2 = [i for i in abs(j2 - nuclearspin):j2 + nuclearspin]
            for f1 in possible_F1, f2 in possible_F2
                # check to see if user defined these levels in level_structure
                hf1 = addhyperfine(copy(e1), f1)
                hf1 in keys(level_structure) || continue
                hf2 = addhyperfine(copy(e2), f2)
                hf2 in keys(level_structure) || continue
                
                # avoid double counting by requiring hf1 to be lower in energy than hf2 
                (level_structure[hf1] < level_structure[hf2]) || continue                
                
                Δf = abs(f1 - f2)
                if multipole == "E1"
                    # ΔF(0↔0)forbidden
                    ((f1 == 0) && (f2 == 0)) && continue
                    # Δ𝐹 and Δ𝐽=0,±1
                    (Δf > 1) && continue
                    # the following should be enforced by users choice of transitions
                    # ΔL(0↔0)forbidden
                    # ((l1 == 0) && (l2 == 0)) && continue
                    # # ΔL = 0, ±1
                    # (ΔL > 1) && continue
                elseif multipole == "M1"
                    # ΔF(0↔0)forbidden
                    ((f1 == 0) && (f2 == 0)) && continue
                    # Δ𝐹 and Δ𝐽=0,±1
                    (Δf > 1) && continue
                    # the following should be enforced by users choice of transitions
                    # # ΔL = 0
                    # (ΔL != 0) && continue
                elseif multipole == "E2"
                    # (0↔0)forbidden
                    ((f1 == 0) && (f2 == 0)) && continue
                    # (0↔1)forbidden
                    ((f1 == 0) && (f2 == 1)) && continue
                    # (1/2↔1/2)forbidden
                    ((f1 == 1//2) && (f2 == 1//2)) && continue
                    # Δ𝐹 and Δ𝐽=0,±2
                    (Δf > 2) && continue
                    # the following should be enforced by users choice of transitions
                    # # ΔL=±2
                    # (ΔL > 2) && continue
                end
                hf_transitions[(hf1, hf2)] = transition_info
            end
        end
    end
    SpeciesProperties(
            shortname,
            mass,
            charge,
            nuclearspin,
            level_structure,
            hf_transitions,
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

Each element of `selected_sublevels` is a 2-element Tuple `(level::EnergyLevel, sublevels)`, 
with the first element being the name of a level and the second specifying which sublevels 
should be included. Allowed sublevels are those whose magnetic quantum number `m` is in the 
set {`-f`, `-f+1`, `-f+2`, ... `f-1`, `f`}, where `f` is the total angular momentum quantum 
number of `level`. For each `level` specified there are three allowed options to specify the 
set of `sublevels` to include:

* `selected_sublevels::Real`: Includes only one `m = sublevels`
* `selected_sublevels::Vector{Real}`: Includes all sublevels whose magnetic quantum number 
    `m` is in `sublevels`
* `selected_sublevels = "all"`: Includes all allowed sublevels

Omission of a level in `selected_sublevels` will exclude all sublevels.

**Fields**
* `speciesproperties`: A container of species-specific properties
* `selected_sublevels`::Vector{Tuple{String,Real}}: Choose which sublevels of Ion to consider
* `shape`::Vector{Int}: Dimension of the Hilbert space
* `manualshift::Dict`: A dictionary with keys denoting the selected levels and values, 
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
    selected_sublevels::Tuple{Vararg{Tuple{EnergyLevel, Real}}}
    selected_transitions::Dict{  # Dictionary of allowed transitions between energy levels
        Tuple{EnergyLevel, EnergyLevel},
        NamedTuple{(:multipole, :einsteinA), Tuple{String, Float64}}
    }
    shape::Vector{Int}             
    manualshift::Dict{Tuple{EnergyLevel, Union{Rational, Int}}, Real}
    ionnumber::Union{Int, Missing}           # remove
    ionposition::Union{Real, Missing}        # remove

    # constructors (overrides default)
    function IonInstance{Species}(
        properties,
        selected_sublevels=nothing,
        manualshift=Dict()
    ) where {Species<:Any}
        sublevels = _construct_sublevels(selected_sublevels, properties)
        selected_transitions = properties.transitions
        shape = [length(sublevels)]
        manualshift_full = _construct_manualshift(manualshift, sublevels)
        return new{Species}(
            properties,
            Tuple(sublevel for sublevel in sublevels),  # Store in an immutable container
            selected_transitions,
            shape,
            manualshift_full,
            missing,          # remove
            missing           # remove
        )
    end

    function IonInstance{Species}(
            speciesproperties,
            sublevels,
            selected_transitions,
            shape,
            manualshift,
            ionnumber,
            ionposition
        ) where {Species<:Any}
            sublevels = deepcopy(sublevels)
            shape = copy(shape)
            manualshift = deepcopy(manualshift)
        return new{Species}(
            speciesproperties,
            sublevels,
            selected_transitions,
            shape,
            manualshift,
            ionnumber,       # remove
            ionposition      # remove
        )
    end
end

function Base.setproperty!(I::Ion, p::Symbol, x)
    if p == :selected_sublevels
        sublevels = _construct_sublevels(x, I.speciesproperties)
        shape = [length(sublevels)]
        Base.setfield!(I, :selected_sublevels, Tuple(sublevel for sublevel in sublevels))
        Base.setfield!(I, :shape, shape)
    elseif p == :shape
        error("setfield!: User shouldn't change shape field directly.")
        return
    else 
        error("setfield!: User shouldn't change this field.")
    end
end


#############################################################################################
# Functions for constructing ion structs
#############################################################################################

function _construct_sublevels(selected_sublevels, properties)
    full_level_structure = properties.level_structure

    # If selected_sublevels is blank, use the default selection. If it is "all", use all 
    # sublevels.
    if isnothing(selected_sublevels) || selected_sublevels == "all" 
        selected_sublevels = [(sublevel, "all") for sublevel in keys(full_level_structure)]
    end

    # Construct the list of sublevels
    sublevels = []
    for manifold in selected_sublevels
        # Ensure that the level exists in the ion
        level = manifold[1]
        @assert level in keys(full_level_structure) "invalid level $level"
        # level_structure = full_level_structure[level]

        selectedms = manifold[2]
        # Add chosen sublevels
        if isnothing(level.f) # this means there is no hyperfine structure
            f = level.j
        else
            f = level.f
        end

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
    return sublevels
end

function _construct_manualshift(manualshift, sublevels)
    manualshift_full = Dict()
    for sublevel in sublevels
        manualshift_full[sublevel] =
            (haskey(manualshift, sublevel) ? manualshift[sublevel] : 0.0)
    end
    return manualshift_full
end


#############################################################################################
# Functions to modify ion fields
#############################################################################################

function validatesublevel(I::Ion, sublevel::Tuple{EnergyLevel, Real})
    @assert sublevel in selected_sublevels(I) (
            "Ion does not contain sublevel $sublevel. Use selected_sublevels(Ion) to see \
            list of available sublevels."
        )
    return sublevel[1]  # returns parent energy level
end

"""
    manualshift!(I::Ion, sublevel, shift::Real)

Applies a manual shift `shift` to the chosen `sublevel` of `I` (overwriting any previously \
assigned manual shift).
"""
function manualshift!(I::Ion, sublevel::Tuple{EnergyLevel, Real}, shift::Real)
    validatesublevel(I, sublevel)
    I.manualshift[(sublevel[1], Rational(sublevel[2]))] = convert(Float64, shift)
    return I.manualshift[(sublevel[1], Rational(sublevel[2]))]
end

"""
    manualshift!(I::Ion, manualshift_dict::Dict)

Applies `manualshift(I, sublevel, shift)` to all pairs `sublevel => shift` of the Dict \
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

Returns a vector of all energy levels of `I`.
"""
levels(I::Ion) = collect(keys(I.speciesproperties.level_structure))

"""
    landegj(l::Real, j::Real, s::Real=1//2)

Landé g-factor of an LS-coupled fine structure energy level. The electron g-factor is taken
to be exactly 2.

**args**
* `l`: orbital angular momentum quantum number
* `j`: electron total angular momentum quantum number
* `s`: electronic spin angular momentum quantum number (defaults to 1/2)
"""
function landegj(l::Real, j::Real, s::Real=1//2)
    (3. / 2) + (s * (s + 1) - l * (l + 1)) / (2j * (j + 1))
end

"""
    landegj(E::LS)

Same thing, but obtain quantum numbers from the state [`LS`](@ref) <: [`EnergyLevel`](@ref).
Note: this quantity is only well-defined for LS-coupled states.
"""
function landegj(E::LS)
    @assert isnothing(E.f) "This EnergyLevel has hyperfine structure."
    landegj(E.l, E.j, E.s)
end

"""
    landegf(l::Real, j::Real, f::Real, i::Real, s::Real=1//2)

Landé g-factor of an LS-coupled hyperfine energy level. The electron g-factor is taken
to be exactly 2.

**args**
* `l`: orbital angular momentum quantum number
* `j`: electron total angular momentum quantum number
* `f`: total angular momentum quantum number
* `i`: nuclear spin angular momentum quantum number
* `s`: electronic spin angular momentum quantum number (defaults to 1/2)
"""
function landegf(l::Real, j::Real, f::Real, i::Real, s::Real=1//2)
    landegj(l, j, s) / 2 * (1 + ((j * (j + 1) - i * (i + 1)) / (f * (f + 1))))
end

"""
    landegf(E::LS[, I::Real])

Same thing, but obtain quantum numbers from the state [`LS`](@ref) <: [`EnergyLevel`](@ref).
Note: this quantity is only well-defined for LS-coupled states. If the total nuclear spin
is not specified by `E`, this can be provided with argument `I`.
"""
function landegf(E::LS)
    @assert !isnothing(E.i) "nuclear spin not specified by energy level, please provide explicitly"
    @assert !isnothing(E.f) "no hyperfine structure specified for this energy level"
    landegf(E.l, E.j, E.f, E.i, E.s)
end
landegf(E::LS, I::Real) = landegf(E.l, E.j, E.f, I, E.s)

"""
    landegf(I::Ion, level::EnergyLevel)

`landegf` for `Energylevel` of `I`. Note: this may be different than what is returned by
`landegf(E<:LS[, I::Real])` since the configuration file for `I` might explicitly provide
a value (which will then be the one used in calculations and returned here).
"""
function landegf(I::Ion, level::EnergyLevel)
    properties = speciesproperties(I)
    if !ismissing(properties.gfactors) && haskey(properties.gfactors, level)
        return properties.gfactors[level]
    else
        if isnothing(level.f)
            return landegj(level)
        else
            return landegf(level, nuclearspin(I))
        end
    end
end
landegj(I::Ion, level::EnergyLevel) = landegf(I, level)

"""
    manualshift(I::Ion, sublevel)

Returns the assigned manual shift of `sublevel` of Ion `I`.
"""
function manualshift(I::Ion, sublevel::Tuple{EnergyLevel, Real})
    validatesublevel(I, sublevel)
    return manualshift(I)[sublevel]
end


"""
    zeemanshift(I::Ion, sublevel, B::Real)

Returns the Zeeman shift at a magnetic field of `B` of `sublevel` of `I`.
If `sublevel` has a custom g-factor defined, then this is used. Otherwise, `landegf` or 
`landegj` is used to compute the Landé g-factor. Zeeman frequency shift calculated as 
``Δf = (μ_B/ħ) ⋅ g_f ⋅ B ⋅ m / 2π``
"""
function zeemanshift(I::Ion, sublevel::Tuple{EnergyLevel, Real}, B::Real)
    validatesublevel(I, sublevel)
    properties = speciesproperties(I)
    if !ismissing(properties.nonlinear_zeeman) &&
       haskey(properties.nonlinear_zeeman, sublevel)
        nonlinear = properties.nonlinear_zeeman[sublevel](B)
    else
        nonlinear = 0.0
    end
    if isnothing(sublevel[1].f)
        return zeemanshift(B, landegj(I, sublevel[1]), sublevel[2]) + nonlinear
    else
        return zeemanshift(B, landegf(I, sublevel[1]), sublevel[2]) + nonlinear
    end
end

"""
    zeemanshift(B::Real, g::Real, m::Real)

Returns (μ_B/ħ) ⋅ g ⋅ B ⋅ m / 2π   
"""
zeemanshift(B::Real, g::Real, m::Real) = (μB / ħ) * g * B * m / 2π

zeemanshift(B::Real, l::Real, j::Real, f::Real, m::Real, i::Real, s::Real=1 // 2) =
    zeemanshift(B, landegf(l, j, f, i, s), m)

"""
    energy(I::Ion, sublevel; B=0, ignore_manualshift=false)

Calculates the energy of `sublevel` of `I`. A Zeeman shift can be included by setting a nonzero 
value of `B`. The manual shift can be omitted by setting `ignore_manualshift=true`.
"""
function energy(I::Ion, sublevel::Tuple{EnergyLevel, Real}; B=0, ignore_manualshift=false)
    validatesublevel(I, sublevel)
    E0 = speciesproperties(I).level_structure[sublevel[1]]
    zeeman = zeemanshift(I, sublevel, B)
    manual = (ignore_manualshift ? 0.0 : manualshift(I, sublevel))
    return E0 + zeeman + manual
end

"""
    transitionfrequency(I::Ion, transition::Tuple; B=0, ignore_manualshift=false)

`transition` is a Tuple of two sublevels or levels. Computes the absolute values of the 
difference in energies between `transition[1]` and `transition[2]`. If between sublevels, 
then the Zeeman shift may be included by setting the value of the magnetic field `B`, and 
manual shifts may be omitted by setting `ignore_manualshift=true`.
"""
function transitionfrequency(
        I::Ion, transition::Tuple{T, T}; B=0, ignore_manualshift=false
    ) where {T<:Tuple{EnergyLevel, Union{Rational, Int}}}
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
    availabletransitions(I::Ion)

Returns a list of allowed transitions of `I`.
"""
availabletransitions(I::Ion) = collect(keys(I.selected_transitions))

"""
    subleveltransitions(I::Ion)

Returns all allowed and selected transitions between sublevels of `I` as a vector of 
`Tuple{S, S} where S<:{EnergyLevel, Real}`.
"""
function subleveltransitions(I::Ion)
    transition_list = []
    for transition in availabletransitions(I)
        (E1, E2) = transition
        # this will fail if one level has hyperfine structure and the other does not
        if isnothing(E1.f)
            (f1, f2) = E1.j, E2.j
        else
            (f1, f2) = E1.f, E2.f
        end
        multipole = transitionmultipole(I, transition)
        sublevels1 = [i for i in -f1:f1]
        sublevels2 = [i for i in -f2:f2]
        for m1 in sublevels1, m2 in sublevels2
            Δm = abs(m2 - m1)
            q = parse(Int, multipole[2])
            # No mj = 0 ↔ 0 if ΔJ or ΔF = 0
            (q == 1) && (m1 == 0 && m2 == 0) && (f1 == f2) && continue
            if Δm <= q
                # Only add to list of sublevel transitions if Δm is not larger than the 
                # transition multipole allows (1 for E1, 2 for E2, etc)
                push!(transition_list, ((E1, m1), (E2, m2)))
            end
        end

    end 
    return transition_list
end

"""
    einsteinA(I::Ion, leveltransition::Tuple{<:EnergyLevel, <:EnergyLevel})

Returns Einstein A coefficient corresponding to the transition 
`leveltransition[1] -> leveltransition[2]`. The first level must be the lower level and the 
second must be the upper level.
"""
function einsteinA(I::Ion, leveltransition::Tuple{<:EnergyLevel, <:EnergyLevel})
    transitions = I.speciesproperties.transitions
    @assert leveltransition in keys(transitions) "invalid transition"
    return transitions[leveltransition].einsteinA
end

"""
    transitionmultipole(I::Ion, leveltransition::Tuple{EnergyLevel, EnergyLevel})

Returns the transition multiple (`'E1'`, `'E2'`, etc.) corresponding to the transition 
`leveltransition[1] -> leveltransition[2]`. The first level must be the lower level and the 
second must be the upper level.
"""
function transitionmultipole(I::Ion, leveltransition::Tuple{EnergyLevel, EnergyLevel})
    transitions = I.speciesproperties.transitions
    @assert leveltransition in keys(transitions) "invalid transition"
    return transitions[leveltransition].multipole
end

"""
    lifetime(I::Ion, level::EnergyLevel)

Computes lifetime of `level` by summing the transition rates out of `level`. The sum is taken 
 over all levels that the species may have, rather than the selected levels. I.e. it is 
 computed based on `I.speciesproperties.transitions` instead of `I.selected_transitions`.
"""
function lifetime(I::Ion, level::EnergyLevel)
    sp = I.speciesproperties
    @assert level in keys(sp.level_structure) (
        "Ion species $(typeof(I)) does not contain level $level"
    )
    totaltransitionrate = 0
    for (transition, info) in sp.transitions
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

#############################################################################################
# Matrix Element Calculations
#############################################################################################

"""
    matrixelement(
        ion::Ion, transition::Tuple, I::Real, ϵhat::NamedTuple, 
        khat::NamedTuple, Bhat::NamedTuple=(;z=1)
    )
Computes the matrix elements (units of Hz) between two energy sublevels.

**args**
* `ion`: Ion undergoing transition
* `transition`: Tuple of sublevels between which the transition is being calculated. Must be 
        formatted such that `energy(transition[2]) > energy(transition[1])`
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

function E1matrixelement()
    nothing
end

function M1matrixelement()
    nothing
end

function E2matrixelement()
    nothing
end

function matrixelement(
        ion::Ion,
        transition::Tuple{T, T},
        I::Real,
        ϵhat::NamedTuple,
        khat::NamedTuple,
        Bhat::NamedTuple=(; z=1)
    ) where {T<:Tuple{EnergyLevel, Union{Rational, Int}}}

    SL1 = transition[1]
    SL2 = transition[2]
    L1 = validatesublevel(ion, SL1)
    L2 = validatesublevel(ion, SL2)

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
        println(sublevel)
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

# TODO: Have this print a nicely formatted summary of pertinent speciesproperties
Base.show(io::IO, m::MIME"text/plain", S::SpeciesProperties) = print("$(S.shortname) SpeciesProperties")
Base.print(io::IO, S::SpeciesProperties) = print("$(S.shortname) SpeciesProperties")
