
export validatesublevel,
    set_sublevel_alias!,
    clear_sublevel_alias!,
    clear_all_sublevel_aliases!,
    levels,
    sublevelalias,
    alias2sublevel,
    sublevel2level,
    leveltransitions,
    subleveltransitions

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
