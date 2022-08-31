include("../ion_configurations/_include_ion_configurations.jl")
export ionstate, sigma, ionprojector

"""
    ionstate(object, sublevel)
For `object<:Ion` and `sublevel<:Tuple{String,Real}` (full sublevel name) or `sublevel<:String`
(alias), this returns the ket corresponding to the `Ion` being in the state ``|index⟩``. The
object can also be an `IonConfiguration` or `Trap` instance, in which case ``N`` arguments
should be given in place of `index`, where ``N`` equals the number of ions in the
`IonConfiguration` or `Trap`. This will return the state
``|index₁⟩⊗|index₂⟩⊗...⊗|index\\_N⟩``.

One may also specify `sublevel<:Int`. If `object<:Ion`, this will return the ket given by
``|index⟩ = (0 ... 1 ... 0)ᵀ`` where the nonzero element in the  column vector is located at
`index`.
"""
function ionstate(I::Ion, sublevel::Tuple{String, Real})
    validatesublevel(I, sublevel)
    i = findall(sublevels(I) .== [sublevel])[1]
    return basisstate(I, i)
end
ionstate(I::Ion, sublevelalias::String) = ionstate(I, alias2sublevel(I, sublevelalias))
ionstate(I::Ion, sublevel::Int) = basisstate(I, sublevel)
function ionstate(IC::IonConfiguration, states::Union{Tuple{String, Real}, String, Int}...)
    ions = IC.ions
    L = length(ions)
    @assert L ≡ length(states) "wrong number of states"
    return tensor([ionstate(ions[i], states[i]) for i in 1:L]...)
end
ionstate(T::Trap, states::Union{Tuple{String, Real}, String, Int}...) =
    ionstate(T.configuration, states...)

"""
    sigma(ion::Ion, ψ1::sublevel[, ψ2::sublevel])
Returns ``|ψ1\\rangle\\langle ψ2|``, where ``|ψ_i\\rangle`` corresponds to the state
returned by `ion[ψᵢ]`.

If ψ2 is not given, then ``|ψ1\\rangle\\langle ψ1|`` is returned.
"""
sigma(ion::Ion, ψ1::T, ψ2::T) where {T <: Union{Tuple{String, Real}, String, Int}} =
    sparse(projector(ion[ψ1], dagger(ion[ψ2])))
sigma(ion::Ion, ψ1::Union{Tuple{String, Real}, String, Int}) = sigma(ion, ψ1, ψ1)

"""
    ionprojector(obj, sublevels...; only_ions=false)
If `obj<:IonConfiguration` this will return ``|ψ₁⟩⟨ψ₁|⊗...⊗|ψ\\_N⟩⟨ψ\\_N|⊗𝟙``
where ``|ψᵢ⟩`` = `obj.ions[i][sublevels[i]]` and the identity operator ``𝟙`` is over all of the
COM modes considered in `obj`.

If `only_ions=true`, then the projector is defined only over the ion subspace.

If instead `obj<:Trap`, then this is the same as `obj = Trap.configuration`.
"""
function ionprojector(
    IC::IonConfiguration,
    sublevels::Union{Tuple{String, Real}, String, Int}...;
    only_ions = false
)
    ions = IC.ions
    L = length(ions)
    @assert L ≡ length(sublevels) "wrong number of sublevels"
    modes = get_vibrational_modes(IC)
    observable = tensor([projector(ions[i][sublevels[i]]) for i in 1:L]...)
    if !only_ions
        for mode in modes
            observable = observable ⊗ one(mode)
        end
    end
    return observable
end
function ionprojector(
    T::Trap,
    sublevels::Union{Tuple{String, Real}, String, Int}...;
    only_ions = false
)
    return ionprojector(T.configuration, sublevels..., only_ions = only_ions)
end
