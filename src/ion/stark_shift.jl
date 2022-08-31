export set_stark_shift!, zero_stark_shift!, stark_shift

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

"""
    stark_shift(I::Ion, sublevel)
Returns the assigned stark shift of `sublevel` of Ion `I`.
"""
function stark_shift(I::Ion, sublevel::Tuple{String, Real})
    validatesublevel(I, sublevel)
    return stark_shift(I)[sublevel]
end
stark_shift(I::Ion, alias::String) = stark_shift(I, alias2sublevel(I, alias))
