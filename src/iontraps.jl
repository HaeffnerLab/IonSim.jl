using LinearAlgebra: eigen
using NLsolve: nlsolve
using Statistics: mean
using Optim
using .PhysicalConstants: e, ϵ₀
using Plots

export IonTrap,
    ions,
    linear_equilibrium_positions,
    LinearChain,
    comfrequencies,
    selectedmodes,
    full_normal_mode_description,
    modes,
    xmodes,
    ymodes,
    zmodes,
    modecutoff!,
    length

"""
    IonTrap
A physical configuration of ions. Stores a collection of ions and information about their
relative center of mass motion.
"""
abstract type IonTrap end

# required functions
ions(I::IonTrap) = I.ions

#############################################################################################
# LinearChain - a linear Coulomb crystal
#############################################################################################

"""
    LinearChain(;
            ions::Vector{Ion}, comfrequencies::NamedTuple{(:x,:y,:z)},
            selectedmodes::NamedTuple{(:x,:y,:z), N::Int=10
        )

Generates and stores all of the information necessary to describe a collection of ions trapped 
in a 3D harmonic potential and forming a linear coulomb crystal. Relevant to calculating the
normal mode structure of the linear chain is the charge, mass  and ordering of the ions.

**user-defined fields**
* `ions::Vector{Ion}`: A list of ions that compose the linear Coulomb crystal
* `comfrequencies::NamedTuple{(:x,:y,:z)`: Describes the COM frequencies 
    `(x=ν_x, y=ν_y, z=ν_z)`. The ``z``-axis is taken to be along the crystal's axis of 
    symmetry. Note that these only correspond to true COM modes if all ions have the same
    mass, otherwise this is a misnomer and corresponds to the smallest (largest) eigenfrequency
    mode in the axial (radial) directions.
* `selectedmodes::NamedTuple{(:x,:y,:z)}`:  e.g. `selectedmodes=(x=[1], y=[2], z=[1,2])`.
    Specifies the axis and a list of integers which correspond to the ``i^{th}`` farthest
    mode away from the COM (see note above about meaning of "COM") for that axis. For example, 
    `selectedmodes=(z=[2])` would specify the axial stretch mode. These are the modes that will 
    be modeled in the chain.Note: `selectedmodes=(x=[],y=[],z=[1])`, `selectedmodes=(y=[],z=[1])`
    and `selectedmodes=(;z=[1])` are all acceptable and equivalent.
* `N::Int=10`: The Hilbert space dimension for each normal mode.
**derived fields**
* `ionpositions::Vector{<:Real}`: The relative positions of the ions in meters.  
"""
struct LinearChain <: IonTrap  # Note: this is not a mutable struct
    ions::Vector{<:Ion}
    comfrequencies::NamedTuple{(:x, :y, :z)}
    selectedmodes::NamedTuple{(:x, :y, :z), Tuple{Vararg{Vector{VibrationalMode}, 3}}}
    ionpositions::Vector{<:Real}
    N::Int
    function LinearChain(; ions, comfrequencies, selectedmodes, N=10, ionpositions=[])
        num_ions = length(ions)
        selectedmodes = _construct_vibrational_modes(selectedmodes)
        warn = nothing
        for i in 1:num_ions, j in (i+1):num_ions
            if ions[j] ≡ ions[i]
                ions[j] = copy(ions[i])
                if isnothing(warn)
                    warn = "Some ions point to the same thing. Making copies."
                    @warn warn
                end
            end
        end
        normalmodes = full_normal_mode_description(ions, comfrequencies)
        vm = (
            x=Vector{VibrationalMode}(undef, 0),
            y=Vector{VibrationalMode}(undef, 0),
            z=Vector{VibrationalMode}(undef, 0)
        )
        axes = [x̂, ŷ, ẑ]
        for (i, modes) in enumerate(selectedmodes), mode in modes
            push!(vm[i], VibrationalMode(normalmodes[i][mode]..., axis=axes[i], N=N))
        end
        l = linear_equilibrium_positions(num_ions)
        kz = searchfor_trappingpotential_parameters(ions, comfrequencies)[1]
        l0 = (e^2 / (4π * ϵ₀ * kz))^(1/3)
        ionpositions = l * l0
        for (i, ion) in enumerate(ions)
            Core.setproperty!(ion, :ionnumber, i)
            Core.setproperty!(ion, :ionposition, ionpositions[i])
        end
        return new(ions, comfrequencies, vm, ionpositions)
    end
end

#############################################################################################
# LinearChain - helper functions
#############################################################################################

_sparsify!(x, eps) = @. x[abs(x)<eps] = 0

#= 
Computes the equilibrium positions for a linear chain of N ions with individual charges
Q[i] (in units of fundamental charge). This equilibrium occurs when global confining 
pseudopotential balances inter-ion Coulomb repulsion for all ions.
=#
function linear_equilibrium_positions(N::Int, Q::Vector{Int})
    @assert !(1 in (Q .< 0) && 1 in (Q .> 0)) (
        "Can't have negative and positive charges in the same chain."
    )
    function f!(F, x, N)
        for i in 1:N
            F[i] = (
                Q[i] * x[i] 
                - sum([Q[i] * Q[j] / (x[i] - x[j])^2 for j in 1:(i-1)]) 
                + sum([Q[i] * Q[j] / (x[i] - x[j])^2 for j in (i+1):N])
            )
        end
    end

    # equilibrium positions follow an empirical relationship used for picking starting point
    # (see ref)
    if isodd(N)
        initial_x = [(2.018 / N^0.559) * i for i in (-N÷2):(N÷2)]
    else
        initial_x =
            [(2.018 / N^0.559) * i for i in filter(x -> x ≠ 0, collect((-N÷2):(N÷2)))]
    end
    sol = nlsolve((F, x) -> f!(F, x, N), initial_x, method = :newton)
    # clean up values
    sol = round.(sol.zero, digits=6)
    for i in eachindex(sol)  # remove annoying signed zeros
        if sol[i] == 0
            sol[i] = 0
        end
    end
    return sol
end

linear_equilibrium_positions(ions::Vector{<:Ion}) = linear_equilibrium_positions(
        length(ions), [ion.speciesproperties.charge for ion in ions]
    )

linear_equilibrium_positions(N::Int) = linear_equilibrium_positions(N, ones(Int, N))

#=
This is the Jacobian for the force of the ions in the chain due to all of the other ions.
Diagonalize it to find the normal modes of motion.
=#
function diagonalize_Kij(
        M_actual::Vector{<:Real}, Q::Vector{<:Int}, axis::NamedTuple{(:x, :y, :z)}, 
        kz::Real, P::Real, kr::Real; 
        optimize=false
    )
    N = length(M_actual)
    maxM = maximum(M_actual)
    M = M_actual / maxM  # scale M_actual, so masses are order unity
    # P *= maxM  # scale P by Mtot since it gets divided by M[i]

    #= 
    kz = ka - kr/2
    kxⱼ = P/Mⱼ - ka/2 + kr   = -kz/2 + 3kr/4
    kyⱼ = P/Mⱼ - ka/2 - kr/2 = -kz/2 - 3kr/4

    Where P is defined s.t. the pseudopotential frequency for jth ion is (√(Qⱼ * P)/Mⱼ) * kz.
    And ka is defined  s.t. the axial DC confining potential frequency is √(Qⱼ * ka / Mⱼ) * kz.
    (Similar for kr, but this is an additional DC quadrupole that is confining in radial direction).
    =#
    if axis == x̂
        a = 1
        kDC =  -kz/2 + 3kr/4  # kx
    elseif axis == ŷ
        a = 1
        kDC = -kz/2 - 3kr/4 # ky
    elseif axis == ẑ
        a = -2
        P = 0
        kDC = kz  # kz
    end

    l = linear_equilibrium_positions(N, Q)
    K = Array{Real}(undef, N, N)
    for i in 1:N, j in 1:N
        if i ≡ j
            K[i, j] = Q[i] * (P/M[i] + kDC) / kz - a * sum([Q[i] * Q[p] / abs(l[j] - l[p])^3 for p in 1:N if p != j])
        else
            K[i, j] = a * Q[i] * Q[j] / abs(l[j] - l[i])^3
        end
        # We switch to mass-scaled coordinates as is standard for linear spring-mass systems.
        K[i, j] /= √(M[i] * M[j])
    end

    bvalues, bvectors = eigen(K)
    bcleanedvalues = []
    # When arg for sqrt is negative, we still want to return something real so the optimizer will work.
    # But we need to make sure it returns a *definitely incorrect* but not super divergent value.
    # This method seems to work, but seems a bit sketchy.
    for value in bvalues
        if kz > 0 && value > 0
            push!(bcleanedvalues, sqrt.(value * kz / maxM))
            # push!(bcleanedvalues, sqrt.(value * kz))
        else
            # we make it negative so it is *definitely incorrect*
            push!(bcleanedvalues, -sqrt.(abs(value * kz / maxM)))
            # push!(bcleanedvalues, -sqrt.(abs(value * kz)))
        end
    end

    # When optimizing, we are just trying to match the computed eigenfrequencies with the target 
    # eigenfrequencies, so no need to compute the eigenvectors
    if optimize
        if axis == ẑ
            return minimum(bcleanedvalues)
        else
            return maximum(bcleanedvalues)
        end
    end

    for i in 1:N
        v = view(bvectors, :, i)
        _sparsify!(v, 1e-6)
        for j in 1:N
            v ./= M[j]  
        end
        # normalization takes care of any remaining global scalings introduced for convenience
        normalize!(v)
        v .*= sign(v[1])
    end

    sortedresults = sort([(bcleanedvalues[i], bvectors[:, i]) for i in 1:length(bcleanedvalues)])
    if axis == ẑ
        # axial modes are sorted from lowest to highest eigenfrequency
        return sortedresults
    else
        # this ordering is reversed for radial modes
        return reverse(sortedresults)
    end
end

diagonalize_Kij(
        ions::Vector{<:Ion}, axis::NamedTuple{(:x, :y, :z)}, 
        kz::Real, P::Real, kr::Real; 
        optimize=false
    ) = diagonalize_Kij(
            [mass(ion) for ion in ions], [ion.speciesproperties.charge for ion in ions],
            axis, kz, P, kr, optimize=optimize
        )

#=
This performs the same function as diagonalize_Kij, but for a homogeneous chain where
all of the ions have the same mass. In that case, one can write the Jacobian Kij in
terms of the characteristic_frequencies, which, in this case are the uniqueCOM frequencies.
Then we can just diagonalize Kij directly instead of having to optimize.
=#
function diagonalize_Kij_for_homogeneous_chain(
        N::Int, com::NamedTuple{(:x, :y, :z)}, axis::NamedTuple{(:x, :y, :z)}
    )
    axis == ẑ ? a = 2 : a = -1
    l = linear_equilibrium_positions(N)
    axis_dict = Dict([(x̂, :x), (ŷ, :y), (ẑ, :z)])
    sa = axis_dict[axis]
    β = com[sa] / com.z
    A = Array{Real}(undef, N, N)
    for n in 1:N, j in 1:N
        if n ≡ j
            A[n, j] = β^2 + a * sum([1 / abs(l[j] - l[p])^3 for p in 1:N if p != j])
        else
            A[n, j] = -a / abs(l[j] - l[n])^3
        end
    end
    a = eigen(A)
    for i in a.values
        @assert i > 0 """
            ($(axis_dict[axis])=$(com[sa]), z=$(com.z)) outside stability region for $N ions
        """
    end
    a1 = [sqrt(i) * com.z for i in a.values]
    a2 = a.vectors
    for i in 1:size(a2, 2)
        v = view(a2, :, i)
        _sparsify!(v, 1e-6)
        v .*= sign(v[1])
    end
    if axis == ẑ
        return [(a1[i], a2[:, i]) for i in 1:length(a1)]
    else
        a = reverse([(a1[i], a2[:, i]) for i in 1:length(a1)])
        return a
    end
end


#=
Given some list of target_eigenfrequencies obtained from the user and some list of
guesses P, ka, kr for the parameters that define the trapping potential, this will
use P, ka, kr to compute the corresponding eigenfrequencies (definin the normal mode 
structure) and then return a quantity describing how far the computed eigenfrequencies 
are from the targets.
=#
function sequentially_diagonalize(
        M::Vector{<:Real}, Q::Vector{<:Int}, target_eigenfrequencies::NamedTuple{(:x, :y, :z)}, 
        kz::Real, P::Real, kr::Real, just_y=false
    )
    N = length(M)
    axes = just_y ? [ŷ] : [x̂, ŷ]
    axes_dict = Dict(x̂=>:x, ŷ=>:y)
    computed_eigenfrequencies = []
    for i in eachindex(axes)
        axis = axes[i]
        computed_eigenfrequency = diagonalize_Kij(M, Q, axis, kz, P, kr, optimize=true)
        push!(computed_eigenfrequencies, computed_eigenfrequency)
    end
    difference = 0
    for i in eachindex(computed_eigenfrequencies)
        difference += abs(computed_eigenfrequencies[i] - target_eigenfrequencies[axes_dict[axes[i]]])
    end
    return difference
end

sequentially_diagonalize(
        ions::Vector{<:Ion}, target_eigenfrequencies::NamedTuple{(:x, :y, :z)}, 
        kz::Real, P::Real, kr::Real,
        just_y=false
    ) = sequentially_diagonalize(
            [mass(ion) for ion in ions], [ion.speciesproperties.charge for ion in ions],
            target_eigenfrequencies, kz, P, kr, just_y=just_y
        )

function searchfor_trappingpotential_parameters(
        M::Vector{<:Real}, Q::Vector{<:Int}, target_eigenfrequencies::NamedTuple{(:x, :y, :z)}
    )
    maxM = maximum(M)
    value = diagonalize_Kij(M, Q, ẑ, maxM, 0, 0; optimize=true)
    kz = (target_eigenfrequencies.z / value)^2 * maxM
    lowerbounds = [-1, -1]
    upperbounds = [Inf, Inf]
    difference = 0
    if target_eigenfrequencies.x != target_eigenfrequencies.y
        initial = [10kz, kz/2]
        res = optimize(
            x -> sequentially_diagonalize(M, Q, target_eigenfrequencies, kz, x[1], x[2]),
            lowerbounds, upperbounds, 
            initial,
            Fminbox(NelderMead()),
            Optim.Options(g_abstol=1e-12, g_reltol=1e-6),
        )
        xfreqs = [i[1] for i in diagonalize_Kij(M, Q, x̂, kz, res.minimizer...)]
        yfreqs = [i[1] for i in diagonalize_Kij(M, Q, ŷ, kz, res.minimizer...)]
        difference += abs(target_eigenfrequencies.x - xfreqs[1])
        difference += abs(target_eigenfrequencies.y - yfreqs[1])
        res = [kz, res.minimizer...]
    else
        initial = [10kz]
        res = optimize(
            x -> sequentially_diagonalize(M, Q, target_eigenfrequencies, kz, x[1], 0),
            lowerbounds[1:1], upperbounds[1:1], 
            initial[1:1],
            Fminbox(NelderMead()),
            Optim.Options(g_abstol=1e-6, g_reltol=1e-6),
        )
        xfreqs = [i[1] for i in diagonalize_Kij(M, Q, x̂, kz, res.minimizer[1], 0)]
        yfreqs = xfreqs
        difference += abs(target_eigenfrequencies.x - xfreqs[1])
        res = [kz, res.minimizer..., 0]
    end

    only_positive_eigenvalues = all(vcat(xfreqs, yfreqs) .> 0)
    only_nonnegative_trapping_parameters = all(res .>= 0)
    error_string = (
    "Solver failed to find a normal mode structure compatible with `characteristicfrequencies`."
    )
    # make sure optimal values are actually close to target_eigenfrequencies
    @assert difference < 1e-3 error_string * " [1]"
    # make sure there are no negative eigenfrequencies
    @assert only_positive_eigenvalues error_string * " [2]"
    # enfore convention that all trapping potentials are nonnegative
    @assert only_nonnegative_trapping_parameters error_string * " [3]"
    # make sure kz > kr, otherwise there is not confinment in z-direction
    @assert res[1] > res[3] error_string * " [4]"

    return res  # This is kz, P, kr
end
    
searchfor_trappingpotential_parameters(
        ions::Vector{<:Ion}, target_eigenfrequencies::NamedTuple{(:x, :y, :z)}
    ) = searchfor_trappingpotential_parameters(
            [mass(ion) for ion in ions], [ion.speciesproperties.charge for ion in ions],
            target_eigenfrequencies
        )


#############################################################################################
# Object fields
#############################################################################################

"""
    comfrequencies(chain::LinearChain)
Returns `chain.comfrequencies`
"""
comfrequencies(chain::LinearChain) = chain.comfrequencies

"""
    selectedmodes(chain::LinearChain)
Returns `chain.selectedmodes`
"""
selectedmodes(chain::LinearChain) = chain.selectedmodes

"""
    full_normal_mode_description(chain::LinearChain)<:NamedTuple{(:x,:y,:z)}
For each axis, this contains an array of tuples where the first element is a vibrational
frequency [Hz] and the second element is a vector describing the participation of each ion
at that vibrational frequency (i.e. the normal mode eigenvector corresponding to that 
eigenfrequency).
"""
function full_normal_mode_description(ions::Vector{<:Ion}, comfreqs::NamedTuple{(:x, :y, :z)})
    trappingparams = searchfor_trappingpotential_parameters(ions, comfreqs)
    normalmodes = (
        x=diagonalize_Kij(ions, x̂, trappingparams...),
        y=diagonalize_Kij(ions, ŷ, trappingparams...),
        z=diagonalize_Kij(ions, ẑ, trappingparams...)
    )
    return normalmodes
end

full_normal_mode_description(lc::LinearChain) = full_normal_mode_description(
    ions(lc), comfrequencies(lc)
)

"""
    modes(lc::LinearChain)
Returns an array of all of the selected `VibrationalModes` in the `LinearChain`.
The order is `[lc.x..., lc.y..., lc.z...]`.
"""
function modes(lc::LinearChain)
    return collect(Iterators.flatten(lc.selectedmodes))
end

"""
    xmodes(lc::LinearChain)
Returns an array of all of the selected `VibrationalModes` in the x-direction in the 
`LinearChain`.
"""
xmodes(lc::LinearChain) = selectedmodes(lc).x
"""
    ymodes(lc::LinearChain)
Returns an array of all of the selected `VibrationalModes` in the y-direction in the 
`LinearChain`.
"""
ymodes(lc::LinearChain) = selectedmodes(lc).y
"""
    zmodes(lc::LinearChain)
Returns an array of all of the selected `VibrationalModes` in the z-direction in the 
`LinearChain`.
"""
zmodes(lc::LinearChain) = selectedmodes(lc).z

"""
    modecutoff!(lc::LinearChain, N::Int)
Sets the upper bound of the Hilbert space of all `VibrationalMode`s in `lc` to be the Fock 
state `N`.
"""
function modecutoff!(lc::LinearChain, N::Int)
    for mode in modes(lc)
        modecutoff!(mode, N)
    end
end

"""	
    basis(chain::LinearChain)::CompositeBasis
Returns the composite basis describing the Hilbert space for `chain`.

Order is ``ion₁ ⊗ ion₂ ⊗ ... ⊗ ion_N ⊗ mode₁ ⊗ mode₂ ⊗ ... ⊗ mode_N``, where the ion
bases are ordered according to the order in `ions(chain)` and the vibrational modes are
ordered according to the order in `[xmodes(chain), ymodes(chain), zmodes(chain)]`.
"""
function basis(chain::LinearChain)
    return tensor(ions(chain)..., xmodes(chain)..., ymodes(chain)..., zmodes(chain))
end

function Base.print(lc::LinearChain)
    println("$(length(lc.ions)) ions")
    println("com frequencies: $(lc.comfrequencies)")
    return println("selected vibrational modes: $(lc.selectedmodes)")
end

function Base.show(io::IO, lc::LinearChain)  # suppress long output
    return print(io, "LinearChain($(length(lc.ions)) ions)")
end

# Takes e.g. (y=[1]) to (x=[], y=[1], z=[])
function _construct_vibrational_modes(x)
    k = collect(keys(x))
    xyz = [:x, :y, :z]
    @assert isnothing(findfirst(x -> x ∉ xyz, k)) (
        "keys of `selectedmodes` must be `:x`, `:y` or `:z`"
    )
    indxs = findall(x -> x ∉ k, xyz)
    values = []
    for i in 1:3
        if i in indxs
            push!(values, Int[])
        else
            push!(values, x[xyz[i]])
        end
    end
    return (; zip(xyz, values)...)
end

characteristic_length_scale(M::Real, ν::Real) = (e^2 / (4π * ϵ₀ * M * (2π * ν)^2))^(1/3)

Base.length(lc::LinearChain) = length(lc.ions)