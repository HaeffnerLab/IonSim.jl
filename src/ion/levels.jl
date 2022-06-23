
export validatesublevel,
    set_sublevel_alias!,
    clear_sublevel_alias!,
    clear_all_sublevel_aliases!,
    levels,
    sublevelalias,
    alias2sublevel,
    sublevel2level,
    energy,
    quantumnumbers

"""
    levels(I::Ion)
Returns array of all energy levels of `I`.
"""
levels(I::Ion) = unique([sublevel[1] for sublevel in sublevels(I)])

"""
    validatesublevel(I::Ion, sublevel)
Asserts that the sublevel given is in the sublevels of the Ion.
Can pass in the alias or a Tuple{String, Real}.
"""
validatesublevel(I::Ion, sublevel::Tuple{String, Real}) =
    @assert sublevel in sublevels(I) "Ion does not contain sublevel $sublevel. Use sublevels(Ion) to see list of available sublevels."
validatesublevel(I::Ion, alias::String) = validatesublevel(I, alias2sublevel(I, alias))

"""
    set_sublevel_alias!(I::Ion, sublevel::Tuple{String,Real}, alias::String)
Assigns an alias `alias` to `sublevel` of `I`.
This then allows one to pass `alias` in place of `sublevel` (for `I` only) into any function which accepts a sublevel as an argument.
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
