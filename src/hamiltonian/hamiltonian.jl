using SparseArrays: rowvals, nzrange, spzeros, spdiagm, findnz, rowvals, nonzeros
using FunctionWrappers: FunctionWrapper
using PolynomialRoots: roots
using QuantumOptics: SparseOperator, embed

export hamiltonian

"""
    hamiltonian(
            T::Trap; timescale::Real=1e-6, lamb_dicke_order::Union{Vector{Int},Int}=1,
            rwa_cutoff::Real=Inf, displacement="truncated", time_dependent_eta=false
        )
Constructs the Hamiltonian for `T` as a function of time. Return type is a function
`h(t::Real, ψ)` that, itself, returns a `QuantumOptics.SparseOperator`.

**args**
* `timescale`: e.g. a value of 1e-6 will take time to be in ``\\mu s``
* `lamb_dicke_order`: Only consider terms that change the phonon number by up to this value.
    If this is an `Int`, then the cutoff is applied to all modes. If this is a `Vector{Int}`,
    then `lamb_dicke_order[i]` is applied to the iᵗʰ mode, according to the order in
    `T.basis`.
    Note: this isn't quite the same thing as the Lamb-Dicke approximation since setting
    `lamb_dicke_order=1` will retain, for example, terms proportional to ``a^\\dagger a ``.
* `rwa_cutoff`: drop terms in the Hamiltonian that oscillate faster than this cutoff.
* `displacement`: This can be either `"truncated"`(default) or `"analytic"`.

   When an atom is irradiated, both the atom's energy and its momentum will generally be
   affected. For an atom in a harmonic potential, the exchange of momentum can be modeled as
   a displacement operation ``D(α=iηe^{-iνt}) = exp[αa^† - α^*a]``, where ``η`` is the
   Lamb-Dicke parameter, which can be described equivalently as either being proportional to
   the square root of the ratio of the recoil frequency with the ground state energy of the
   atom's motion or as the ratio of the spread of the ground state wavefunction to the
   wavelength of the laser.

   When `"truncated"` is selected, the matrix elements of ``D(α)`` are computed by
   constructing ``α^* a, αa^†`` in a truncated basis (according to the dimension specified in
   your model) and then exponentiating their difference. This has the advantage, amongst
   other things, of guaranting unitarity.

   If `"analytic"` is selected, then the matrix elements are computed assuming an infinite-
   dimensional Hilbert space.

   For small displacements (``η ≪ N``, where ``N`` is the dimension of the motion's Hilbert
   space), both of these methods will be good approximations.
* `time_dependent_eta::Bool`: In addition to impacting the vibrational subspace directly, a
   change in the trap frequency, ``δν``, will also change the Lamb-Dicke parameter. Since
   typically ``δν≪ν``, this effect will be small ``η ≈ η₀(1 + δν/2ν)`` and doesn't warrant
   the additional computational resources needed to calculate and update it in time. In this
   case, we can set `time_dependent_eta=false` (default), which will set ``η(t) = η₀``.

"""
function hamiltonian(
    T::Trap;
    timescale::Real = 1e-6,
    lamb_dicke_order::Union{Vector{Int}, Int} = 1,
    rwa_cutoff::Real = Inf,
    displacement::String = "truncated",
    time_dependent_eta::Bool = false
)
    return hamiltonian(
        T,
        T.configuration,
        timescale,
        lamb_dicke_order,
        rwa_cutoff,
        displacement,
        time_dependent_eta
    )
end

#############################################################################################
# Hamiltonian for a linear configuration of ions
#############################################################################################

# At each time step, this function updates in-place the 2D array describing the full system
# Hamiltonian.
function hamiltonian(
    T::Trap,
    configuration::LinearChain,
    timescale::Real,
    lamb_dicke_order::Union{Vector{Int}, Int},
    rwa_cutoff::Real,
    displacement::String,
    time_dependent_eta::Bool
)
    b, indxs, cindxs = _setup_base_hamiltonian(
        T,
        timescale,
        lamb_dicke_order,
        rwa_cutoff,
        displacement,
        time_dependent_eta
    )
    aui, gbi, gbs, bfunc, δνi, δνfuncs = _setup_fluctuation_hamiltonian(T, timescale)
    S = SparseOperator(get_basis(T))
    function f(t, ψ)  # a two argument function is required in the QuantumOptics solvers
        @inbounds begin
            @simd for i in 1:length(indxs)
                bt_i, conj_bt_i = b[i](t)::Tuple{ComplexF64, ComplexF64}
                @simd for j in 1:length(indxs[i])
                    i1, i2 = indxs[i][j]
                    S.data[i1, i2] = bt_i
                    if i1 != i2
                        S.data[i2, i1] = conj(bt_i)
                        if length(cindxs[i]) != 0
                            flag = cindxs[i][1][1]
                            i3, i4 = cindxs[i][j + 1]
                            if flag == -1
                                S.data[i3, i4] = -conj_bt_i
                                S.data[i4, i3] = -conj(conj_bt_i)
                            else
                                S.data[i3, i4] = conj_bt_i
                                S.data[i4, i3] = conj(conj_bt_i)
                            end
                        end
                    end
                end
            end
            if length(gbi) == 0 && length(δνi) == 0
                return S
            else
                @simd for indx in aui
                    S.data[indx, indx] = complex(0.0)
                end
                @simd for i in 1:length(gbi)
                    zeeman_t = bfunc(t)::Float64
                    @simd for j in 1:length(gbi[i])
                        indx = gbi[i][j]
                        S.data[indx, indx] += zeeman_t * gbs[i][j]
                    end
                end
                @simd for i in 1:length(δνi)
                    δν_t = δνfuncs[i](t)::Float64
                    @simd for n in 1:length(δνi[i])
                        @simd for indx in δνi[i][n]
                            S.data[indx, indx] += n * δν_t
                        end
                    end
                end
            end
        end
        return S
    end
    return f
end


