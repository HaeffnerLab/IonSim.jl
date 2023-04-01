export TermSymbol,
    LSTerm,
    J₁KTerm,
    @jk_str,
    @ls_str


"""
    TermSymbol

Base type for representing a spectroscopic term.
"""
abstract type TermSymbol end

"""
    LSTerm <: TermSymbol

Stores the same information as a spin-orbit coupled spectroscopic term. 

**args**
`n::Int`: principal quantum number
`l::Int`: orbital quantum number
`s::Union{Int, Rational{Int}}`: total spin angular momentum
`j::Union{Int, Rational{Int}}`: total angular momentum quantum number
"""
@with_kw struct LSTerm <: TermSymbol @deftype Union{Int, Rational{Int}}
    n::Int
    l::Int
    s
    j
end

"""
    J₁KTerm <: TermSymbol    

Stores the same information as an intermediate, J₁K coupled spectroscopic term. 

**args**
`k:::Union{Int, Rational{Int}}`: quantum number describing coupling between J₁ and K
`s:::Union{Int, Rational{Int}}`: total spin quantum number for outermost valence electrons
`j::Union{Int, Rational{Int}}`: total electronic angular momentum (spin and orbital)
"""
@with_kw struct J₁KTerm <: TermSymbol @deftype Union{Int, Rational{Int}}
    # It seems that we can't specify a particular n?
    # n
    k
    s
    j
end

"""
    jk_str(T::AbstractString)

Macro for constructing `J₁KTerm`` from string. Format is e.g. ³[3/2]_1/2 and can be 
constructed with `jk"³[3/2]_1/2"`.
"""
macro jk_str(T::AbstractString)
    lb = findfirst('[', T)  # left bracket location
    rb = findfirst(']', T)  # right bracket location
    ul = findfirst('_', T)  # underscore location
    spin_dict = Dict('¹'=>1, '²'=>2, '³'=>3)
    
    error_msg(i) = (
        """
        Problem around value '$(T[i])' in position $i.
        Syntax for J₁K-coupling term is: 
        1st character               : 2s + 1
        characters between brackets : k
        characters after _          : j
        E.g. lowest triplet 3/2 state in ¹⁷¹Yb⁺ is ³[3/2]_1/2
        """
    )
    
    isnothing(ul) && error(error_msg(6))  # no underscore
    S = occursin(T[1], "¹²³")     ? spin_dict[T[1]]     : error(error_msg(1))
    s = iseven(S) ? (S - 1) // 2 : Int((S - 1) / 2)
    
    fraction = T[lb+1:rb-1]
    fraction_tuple = split(fraction, "/")
    k = _checkfraction(fraction_tuple)
    isnothing(k) && error(error_msg(3))
    
    fraction = T[ul+1:end]
    fraction_tuple = split(fraction, "/")
    j = _checkfraction(fraction_tuple)
    isnothing(j) && error(error_msg(7))
    
    return J₁KTerm(k=k, s=s, j=j)
end

"""
    ls_str(T::AbstractString)

Macro for constructing `LSTerm`` from string. Format is e.g. 4²S_1/2 and can be constructed
with `ls"4²S_1/2"`.
"""
macro ls_str(T::AbstractString)
    ul = findfirst('_', T)  # underscore location
    spin_dict = Dict('¹'=>1, '²'=>2, '³'=>3)
    oam_dict = Dict('S'=>0, 'P'=>1, 'D'=>2, 'F'=>3)
    
    error_msg(i) = (
        "
        Problem around value '$(T[i])' in position $i.
        Syntax for LS-coupling term is: 
        1st character       : n
        2nd character       : 2s + 1
        3rd character.      : l
        characters after _  : j
        E.g. ⁴⁰Ca⁺ ground state is 4²S_1/2.
        "
    )
    
    isnothing(ul) && error(error_msg(4))  # no underscore
    n = isdigit(T[1])             ? parse(Int, T[1])    : error(error_msg(1))
    S = occursin(T[2], "¹²³")     ? spin_dict[T[2]]     : error(error_msg(2))
    s = iseven(S) ? (S - 1) // 2 : Int((S - 1) / 2)
    l = occursin(T[ul-1], "SPDF") ? oam_dict[T[ul-1]]   : error(error_msg(3))
    
    fraction = T[ul+1:end]
    fraction_tuple = split(fraction, "/")
    j = _checkfraction(fraction_tuple)
    isnothing(j) && error(error_msg(4))
    
    return LSTerm(n=n, s=s, l=l, j=j)
end

function _checkfraction(fraction)
    if length(fraction) == 1
        return parse(Int, fraction[1])
    elseif length(fraction) == 2
        num, denom = fraction
        isdigit(num[1]) || error(error_msg(4))   
        isdigit(denom[1]) || error(error_msg(4))
        num, denom = parse.(Int, (num, denom))
        return num // denom
    else
        return nothing
    end
end

function Base.print(T::TermSymbol)
    _tostr(x) = typeof(x)<:Int ? "$x" : "$(x.num)/$(x.den)"
    str = "|"
    for sym in fieldnames(typeof(T))
        sym_str = String(sym)
        x = _tostr(getfield(T, sym))
        str *= "$(sym_str)=$x, "
    end
    str = str[1:end-2] * "⟩"
    print(str)
end

function Base.println(T::TermSymbol)
    print(T)
    println()
end
