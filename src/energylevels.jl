export EnergyLevel,
    LS,
    J₁K,
    @jk_str,
    @ls_str,
    addhyperfine


"""
    EnergyLevel

Base type for representing an electronic energy level corresponding to a spectroscopic term.
"""
abstract type EnergyLevel end

"""
    LS <: EnergyLevel

Stores description of an LS-coupled energy level. This is only intended to be instantiated
by the string macro [`@ls_str`](@ref)

**args**
* `n`: principal quantum number
* `l`: electronic orbital angular momentum 
* `s`: electronic spin angular momentum
* `j`: total electronic angular momentum
**optional**
* `f`: total angular momentum (including nuclear spin)
* `str_repr`: a string representation of the level
"""
@with_kw struct LS <: EnergyLevel
    n::Int 
    l::Int; @assert l < n
    s::Union{Int, Rational{Int}}
    j::Union{Int, Rational{Int}}; @assert abs(l-s) <= j <= l+s
    f::Union{Int, Rational{Int}, Nothing}=nothing; # we don't require I, so can't restrict values
    str_repr::Union{String, Nothing}=nothing
end

"""
    J₁K <: EnergyLevel   

Stores description of an intermediate, J₁K coupled energy level. This is only intended to be 
instantiated by the string macro [`@jk_str`](@ref).

**args**
* `k`: quantum number describing angular momentum coupling between J₁ and K
* `s`: total spin quantum number for outermost valence electrons
* `j`: total electronic angular momentum (spin and orbital)
**optional**
* `f`: total angular momentum (including nuclear spin) 
* `str_repr`: a string representation of the level
"""
@with_kw struct J₁K <: EnergyLevel @deftype Union{Int, Rational{Int}}
    # It seems that we can't specify a particular n?
    # n
    k
    s
    j
    f::Union{Int, Rational{Int}, Nothing}=nothing # we don't require I, so can't restrict values
    str_repr::Union{String, Nothing}=nothing
end

"""
    jk_str(T::AbstractString)

Macro for constructing [`J₁K`](@ref) `<: EnergyLevel` from string. 
    
```julia-repl
julia> jk"³[3/2]_1/2"
|k=3/2, s=1, j=1/2⟩
``` 
\n
```julia-repl
julia> jk"³[3/2]_1/2(F=0)"
|k=3/2, s=1, j=1/2, f=0⟩
``` 
"""
macro jk_str(T::AbstractString)
    lb = findfirst('[', T)  # left bracket location
    rb = findfirst(']', T)  # right bracket location
    ul = findfirst('_', T)  # underscore location
    pl = findfirst('(', T)  # parantheses location
    spin_dict = Dict('¹'=>1, '²'=>2, '³'=>3)
    error_msg = "Bad energy level syntax. Proper example: ³[3/2]_1/2"
    
    isnothing(ul) && error(error_msg(6))  # no underscore
    S = occursin(T[1], "¹²³")     ? spin_dict[T[1]]     : error(error_msg)
    s = iseven(S) ? (S - 1) // 2 : Int((S - 1) / 2)
        
    k = _checkfraction(T[lb+1:rb-1])
    isnothing(k) && error(error_msg)
    
    slice = isnothing(pl) ? T[ul+1:end] : T[ul+1:pl-1]
    j = _checkfraction(slice)
    isnothing(j) && error(error_msg)

    if isnothing(pl)
        f = nothing
    else
        f = _checkforhyperfine(T)
    end
    
    return J₁K(k=k, s=s, j=j, f=f, str_repr=T)
end

"""
    ls_str(T::AbstractString)

Macro for constructing `LS <: EnergyLevel` from string. 

```julia-repl
julia> ls"4²S_1/2"
|n=4, l=0, s=1/2, j=1/2⟩
```

```julia-repl
julia> ls"4²S_1/2(F=0)"
|n=4, l=0, s=1/2, j=1/2, f=0⟩
```
"""
macro ls_str(T::AbstractString)
    ul = findfirst('_', T)  # underscore location
    pl = findfirst('(', T)  # parantheses location
    spin_dict = Dict('¹'=>1, '²'=>2, '³'=>3)
    oam_dict = Dict('S'=>0, 'P'=>1, 'D'=>2, 'F'=>3)
    error_msg = "Bad energy level syntax. Proper example: ⁴⁰Ca⁺ ground state is 4²S_1/2."
    
    isnothing(ul) && error(error_msg(4))  # no underscore
    n = isdigit(T[1])             ? parse(Int, T[1])    : error(error_msg)
    S = occursin(T[2], "¹²³")     ? spin_dict[T[2]]     : error(error_msg)
    l = occursin(T[ul-1], "SPDF") ? oam_dict[T[ul-1]]   : error(error_msg)
    s = iseven(S) ? (S - 1) // 2 : Int((S - 1) / 2)
    
    slice = isnothing(pl) ? T[ul+1:end] : T[ul+1:pl-1]
    j = _checkfraction(slice)
    isnothing(j) && error(error_msg)
    
    if isnothing(pl)
        f = nothing
    else
        f = _checkforhyperfine(T)
    end

    return LS(n=n, s=s, l=l, j=j, f=f, str_repr=T)
end

function _checkfraction(fraction)
    fraction = split(fraction, "/")
    if length(fraction) == 1
        return parse(Int, fraction[1])
    elseif length(fraction) == 2
        num, denom = fraction
        isdigit(num[1]) || error(error_msg)   
        isdigit(denom[1]) || error(error_msg)
        num, denom = parse.(Int, (num, denom))
        return num // denom
    else
        return nothing
    end
end

function _checkforhyperfine(s)
    if occursin('(', s) && occursin(')', s)
        n1 = findfirst('(', s)
        n2 = findfirst(')', s)
        ss = split(s[n1+1:n2-1], '=')
        f = length(ss) == 2 ? ss[2] : nothing
    else
        f = nothing
    end
    if !isnothing(f)
        f = _checkfraction(f)
    end
    return f
end

function addhyperfine(T::EnergyLevel, F::Union{Rational, Int})
    @assert isnothing(T.f) "!isnothing(EnergyLevel.f)" 
    try
        F = convert(Int, F)
    catch InexactError
        F = "$(F.num)/$(F.den)"
    end
    return typeof(T)(T.n, T.l, T.s, T.j, F, T.str_repr * "(F=$F)")
end

#TODO: Custom type for sublevels (currently they are stored as Tuple{EnergyLevel, Real})


#############################################################################################
# Overrides of Base functions
#############################################################################################

# we don't care about str_repr
function Base.:(==)(E1::EnergyLevel, E2::EnergyLevel)
    (E1.n == E2.n) && (E1.l == E2.l) && (E1.s == E2.s) && (E1.j == E2.j) && (E1.f == E2.f)
end

function Base.show(io::IO, ::MIME"text/plain", T::EnergyLevel)
    _tostr(x) = typeof(x)<:Int ? "$x" : "$(x.num)/$(x.den)"
    str = isnothing(T.str_repr) ? "|" : "$(T.str_repr): |"
    for sym in fieldnames(typeof(T))
        sym == :str_repr && continue
        sym_str = String(sym)
        x = getfield(T, sym)
        isnothing(x) && continue
        x = _tostr(x)
        str *= "$(sym_str)=$x, "
    end
    str = str[1:end-2] * "⟩"
    println(io, str)
end

function Base.dump(io::IO, T::EnergyLevel; maxdepth=8, indent="")
    if isnothing(T.str_repr)
        print(io, T)
    else
        print(io, T.str_repr)
    end
end
    
function Base.show(io::IO, x::Tuple{Vararg{Tuple{<:EnergyLevel, <:Any}}})
    for t in x
        println(io, t)
    end
end

function Base.show(io::IO, t::Tuple{EnergyLevel, Union{Rational, Int}})
    t1 = t[1]
    t2 = t[2]
    try
        t2 = convert(Int, t2)
    catch InexactError
        t2 = "$(t2.num)/$(t2.den)"
    end
    print(io, "(", t1, ", ", t2, ")")
end

function Base.show(io::IO, ::MIME"text/plain", x::AbstractDict{<:EnergyLevel, <:Any})
    println(io, "OrderedDict{EnergyLevel, Real} with $(length(x)) entries:")
    for pair in x
        println(io, pair)
    end
end



