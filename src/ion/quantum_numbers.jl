
export quantumnumbers

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
