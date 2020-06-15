using QuantumOptics: projector, tensor, SparseOperator, DenseOperator, basisstate, Ket
using LinearAlgebra: diagm
import QuantumOptics: displace, thermalstate, coherentthermalstate, fockstate


export create, destroy, number, displace, coherentstate, coherentthermalstate, fockstate,
       thermalstate, sigma, ionprojector, ionstate


#############################################################################################
# VibrationalMode operators
#############################################################################################

"""
    create(v::VibrationalMode)
returns the creation operator for `v` such that: `create(v) * v[i] = √(i+1) * v[i+1]`.
"""
create(v::VibrationalMode) = SparseOperator(v, diagm(-1 => sqrt.(1:v.N)))

"""
    destroy(v::VibrationalMode)
Returns the destruction operator for `v` such that: `destroy(v) * v[i] = √i * v[i-1]`.
"""
destroy(v::VibrationalMode) = create(v)'

"""
    number(v::VibrationalMode)
Returns the number operator for `v` such that:  `number(v) * v[i] = i * v[i]`.
"""
number(v::VibrationalMode) = SparseOperator(v, diagm(0 => 0:v.N))

"""
    displace(v::VibrationalMode, α::Number)
Returns the displacement operator ``D(α)`` corresponding to `v`.
"""
function displace(v::VibrationalMode, α::Number)
    D = zeros(ComplexF64, v.N+1, v.N+1)
    @inbounds for n in 1:v.N+1, m in 1:v.N+1
        D[n, m] = _Dnm(α, n, m)
    end
    DenseOperator(v, D)
end

"""
    thermalstate(v::VibrationalMode, n̄::Real)
Returns a thermal density matrix with ``⟨a^†a⟩ ≈ n̄``. Note: approximate because we are 
dealing with a finite dimensional Hilbert space that must be normalized.
"""
function thermalstate(v::VibrationalMode, n̄::Real)
    if n̄ == 0
        return v[0] ⊗ v[0]'
    end
    T = 1 / (log((1 / n̄) + 1))
    H = create(v) * destroy(v)
    thermalstate(H, T)
end

"""
    coherentstate(v::VibrationalMode, α::Number)
Returns a coherent state on `v` with complex amplitude ``α``.
"""
function coherentstate(v::VibrationalMode, α::Number)
    # this implementation is the same as in QuantumOptics.jl, but there the function is 
    # restricted to v::FockBasis, so we must reimplement here
    k = zeros(ComplexF64, v.N+1)
    k[1] = exp(-abs2(α) / 2)
    @inbounds for n=1:v.N
        k[n+1] = k[n] * α / √n
    end
    Ket(v, k)
end

"""
    coherentthermalstate(v::VibrationalMode, n̄::Real, α::Number)
Returns a displaced thermal state for `v`. The mean occupation of the thermal state is `n̄` 
, and `α` is the complex amplitude of the displacement.
"""
function coherentthermalstate(v::VibrationalMode, n̄::Real, α::Number)
    d = displace(v, α)
    d * thermalstate(v, n̄) * d'
end

"""
    fockstate(v::VibrationalMode, N::Int)
Returns the fockstate ``|N⟩`` on `v`. 
"""
fockstate(v::VibrationalMode, N::Int) = v[N]

#############################################################################################
# Ion operators
#############################################################################################

"""
    ionstate(object, index)
For `object<:Ion` and `index<:String`, this returns the ket corresponding to the `Ion` being 
in the state ``|index⟩``. The object can also be an `IonConfiguration` or `Trap` instance, in
which case ``N`` arguments should be given in place of `index`, where ``N`` equals the number
of ions in the `IonConfiguration` or `Trap`. This will return the state 
``|index₁⟩⊗|index₂⟩⊗...⊗|index\\_N⟩``.

Rather than `index<:String`, one may also specify `index<:Int`. If `object<:Ion`, this will
return the ket given by ``|index⟩ = (0 ... 1 ... 0)ᵀ`` where the nonzero element in the 
column vector is located at `index`.
"""
function ionstate(I::Ion, state::String)
    s = I.selected_level_structure.keys
    @assert state in s "index not in selected_level_structure: $s"
    i = findall(s .≡ state)[1]
    basisstate(I, i)
end

function ionstate(IC::IonConfiguration, states::Union{String,Int}...)
    ions = IC.ions
    L = length(ions)
    @assert L ≡ length(states) "wrong number of states"
    tensor([ionstate(ions[i], states[i]) for i in 1:L]...)
end

ionstate(T::Trap, states::Union{String,Int}...) = ionstate(T.configuration, states...)
ionstate(I::Ion, level::Int) = basisstate(I, level)

"""
    sigma(ion::Ion, ψ1::Union{String,Int}, ψ2::Union{String,Int})
Returns ``|ψ1\\rangle\\langle ψ2|``, where ``|ψ_i\\rangle`` corresponds to the state
returned by `ion[ψᵢ]`.
"""
sigma(ion::Ion, ψ1::T, ψ2::T) where {T<:Union{String,Int}} = projector(ion[ψ1], dagger(ion[ψ2]))

"""
    ionprojector(obj, states::Union{String,Int}...; only_ions=false)

If `obj<:IonConfiguration` this will return ``|ψ₁⟩⟨ψ₁|⊗...⊗|ψ\\_N⟩⟨ψ\\_N|⊗𝟙`` 
where ``|ψᵢ⟩`` = `obj.ions[i][states[i]]` and the identity operator ``𝟙`` is over all of the 
COM modes considered in `obj`.

If `only_ions=true`, then the projector is defined only over the ion subspace.

If instead `obj<:Trap`, then this is the same as `obj = Trap.configuration`.
"""
function ionprojector(IC::IonConfiguration, states::Union{String,Int}...; only_ions=false)
    ions = IC.ions
    L = length(ions)
    @assert L ≡ length(states) "wrong number of states"
    modes = get_vibrational_modes(IC)
    observable = tensor([projector(ions[i][states[i]]) for i in 1:L]...)
    if !only_ions
        for mode in modes
            observable = observable ⊗ one(mode)
        end
    end
    observable
end

function ionprojector(T::Trap, states::Union{String,Int}...; only_ions=false)
    ionprojector(T.configuration, states..., only_ions=only_ions)
end 