using LinearAlgebra: eigen
using NLsolve: nlsolve
using Statistics: mean
using Optim
using .PhysicalConstants: e, ϵ₀
using Plots
import YAML

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
    length,
    ionpositions,
    LinearChain_fromyaml,
    visualize

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
            ions, comfrequencies::NamedTuple{(:x,:y,:z)},
            selectedmodes::NamedTuple{(:x,:y,:z), N::Int=10
        )

Generates and stores all of the information necessary to describe a collection of ions 
trapped in a 3D harmonic potential and forming a linear coulomb crystal. Relevant to 
calculating the normal mode structure of the linear chain is the charge, mass  and ordering 
of the ions.

**user-defined fields**
* `ions`: An iterable list of ions that compose the linear Coulomb crystal
* `comfrequencies::NamedTuple{(:x,:y,:z)`: Describes the COM frequencies 
    `(x=ν_x, y=ν_y, z=ν_z)`. The z-axis is taken to be along the crystal's axis of 
    symmetry. Note that these only correspond to true COM modes if all ions have the same
    mass, otherwise this is a misnomer and corresponds to the smallest (largest) 
    eigenfrequency mode in the axial (radial) directions. Note: we assume that 
    `ν_x > ν_y > ν_z`.
* `selectedmodes::NamedTuple{(:x,:y,:z)}`:  e.g. `selectedmodes=(x=[1], y=[1, 2:3], z=[:])`.
    Specifies the axis and a list of integers which correspond to the ``i^{th}`` farthest
    mode away from the COM (see note above about meaning of "COM") for that axis. For example, 
    `selectedmodes=(z=[2])` would specify the axial stretch mode. These are the modes that 
    will be modeled in the chain.Note: `selectedmodes=(x=[],y=[],z=[1])`, 
    `selectedmodes=(y=[],z=[1])` and `selectedmodes=(;z=[1])` are all acceptable and 
    equivalent.
* `N::Int=10`: Optionally specify the Hilbert space dimension for each normal mode.
**derived fields**
* `ionpositions::Vector{<:Real}`: The relative positions of the ions in meters.  
"""
struct LinearChain <: IonTrap
    ions::Tuple{Vararg{Ion}}
    comfrequencies::NamedTuple{(:x, :y, :z)}
    selectedmodes::NamedTuple{(:x, :y, :z), Tuple{Vararg{Vector{VibrationalMode}, 3}}}
    ionpositions::Vector{<:Real}

    function LinearChain(;
        ions,
        comfrequencies::NamedTuple,
        selectedmodes::NamedTuple,
        N::Int=10
    )
        try
            ions = collect(ions)
        catch
            AssertionError("`ions` must be an iterable collection of `<:Ion`.") |> throw
        end
        num_ions = length(ions)
        selectedmodes = _construct_vibrational_modes(selectedmodes, num_ions)

        # shouldn't need this after 
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
        l0 = (e^2 / (4π * ϵ₀ * kz))^(1 / 3)
        ionpositions = l * l0
        _setionproperties(ions, ionpositions)
        return new(ions |> Tuple, comfrequencies, vm, ionpositions)
    end

    function LinearChain(ions, vm, ionpositions)
        _setionproperties(ions, ionpositions)
        new(ions |> Tuple, (x=NaN, y=NaN, z=NaN), vm, ionpositions)
    end
end

"""
    LinearChain_fromyaml(ions::Vector{<:Ion}, yaml::String; N::Int=10)

Load normal mode structure from the specified `yaml` file. It is up to the user to enforce 
the physicality of this structure. The yaml file should have the following structure:

````
---
x:
  - frequency: 1e6
    mode: [0.1, 0.5, 0.3, 0.8]
  - frequency: 2e6
    mode: [0.3, 0.6, 0.5, 3]
y:
  - frequency: 8e6
    mode: [1, 1, 1, 1]
ionpositions: [-1, -0.5, -0.25, 5]
````
"""
function LinearChain_fromyaml(; ions::Vector{<:Ion}, yaml::String, N::Int=10)
    num_ions = length(ions)
    axes = ["x", "y", "z"]
    axes_dict = Dict("x" => x̂, "y" => ŷ, "z" => ẑ)
    # comfrequencies = (x=NaN, y=NaN, z=NaN)
    yaml = YAML.load_file(yaml)
    k = keys(yaml)
    @assert "ionpositions" in k ("yml file must have a key 'ionpositions'")
    ionpositions = yaml["ionpositions"]
    @assert length(ionpositions) == num_ions (
        "length(ionpositions) must equal length(ions)"
    )
    vm = Dict("x" => [], "y" => [], "z" => [])
    for axis in axes
        if axis in k
            axismodes = yaml[axis]
            for mode in axismodes
                if !("frequency" in keys(mode))
                    @warn "no frequency specified for a mode in $axis"
                    continue
                elseif !("mode" in keys(mode))
                    @warn "no mode specified for frequency $(mode["frequency"])"
                    continue
                elseif !(length(mode["mode"]) == num_ions)
                    @warn (
                        "length(mode) ≠ length(ions) for $(mode["frequency"]), \
                        $(mode["mode"])"
                    )
                    continue
                end
                push!(vm[axis], (mode["frequency"], mode["mode"]))
            end
        end
    end
    modes = (x=vm["x"], y=vm["y"], z=vm["z"])
    vm = (
        x=Vector{VibrationalMode}(undef, 0),
        y=Vector{VibrationalMode}(undef, 0),
        z=Vector{VibrationalMode}(undef, 0)
    )
    for (i, axis) in enumerate(axes), mode in modes[Symbol(axis)]
        push!(vm[i], VibrationalMode(mode..., axis=axes_dict[axis], N=N))
    end
    println(vm)
    LinearChain(ions, vm, ionpositions)
    # LinearChain(
    #     ions=ions, comfrequencies=comfrequencies, selectedmodes=modes, 
    #     ionpositions=ionpositions
    # )
end

function _setionproperties(ions, ionpositions)
    for (i, ion) in enumerate(ions)
        Core.setproperty!(ion, :ionnumber, i)
        Core.setproperty!(ion, :ionposition, ionpositions[i])
    end
end
#############################################################################################
# LinearChain - helper functions
#############################################################################################

_sparsify!(x, eps) = @. x[abs(x)<eps] = 0

# Computes the equilibrium positions for a linear chain of N ions with individual charges
# Q[i] (in units of fundamental charge). This equilibrium occurs when global confining 
# pseudopotential balances inter-ion Coulomb repulsion for all ions.
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
    sol = nlsolve((F, x) -> f!(F, x, N), initial_x, method=:newton)
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
    length(ions),
    [ion.speciesproperties.charge for ion in ions]
)

linear_equilibrium_positions(N::Int) = linear_equilibrium_positions(N, ones(Int, N))

# This is the Jacobian for the force of the ions in the chain due to all of the other ions.
# Diagonalize it to find the normal modes of motion.
function diagonalize_Kij(
    M_actual::Vector{<:Real},
    Q::Vector{<:Int},
    axis::NamedTuple{(:x, :y, :z)},
    kz::Real,
    P::Real,
    kr::Real;
    optimize=false
)
    N = length(M_actual)
    maxM = maximum(M_actual)
    M = M_actual / maxM  # scale M_actual, so masses are order unity
    # P *= maxM  # scale P by Mtot since it gets divided by M[i]


    # kz = ka - kr/2
    # kxⱼ = P/Mⱼ - ka/2 + kr   = -kz/2 + 3kr/4
    # kyⱼ = P/Mⱼ - ka/2 - kr/2 = -kz/2 - 3kr/4

    # Where P is defined s.t. the pseudopotential frequency for jth ion is (√(Qⱼ * P)/Mⱼ) * kz.
    # And ka is defined  s.t. the axial DC confining potential frequency is 
    # √(Qⱼ * ka / Mⱼ) * kz. (Similar for kr, but this is an additional DC quadrupole that is 
    # confining in radial direction).
    if axis == x̂
        a = 1
        kDC = -kz / 2 + 3kr / 4  # kx
    elseif axis == ŷ
        a = 1
        kDC = -kz / 2 - 3kr / 4 # ky
    elseif axis == ẑ
        a = -2
        P = 0
        kDC = kz  # kz
    end

    l = linear_equilibrium_positions(N, Q)
    K = Array{Real}(undef, N, N)
    for i in 1:N, j in 1:N
        if i ≡ j
            K[i, j] =
                Q[i] * (P / M[i] + kDC) / kz -
                a * sum([Q[i] * Q[p] / abs(l[j] - l[p])^3 for p in 1:N if p != j])
        else
            K[i, j] = a * Q[i] * Q[j] / abs(l[j] - l[i])^3
        end
        # We switch to mass-scaled coordinates as is standard for linear spring-mass systems.
        K[i, j] /= √(M[i] * M[j])
    end

    bvalues, bvectors = eigen(K)
    bcleanedvalues = []
    # When arg for sqrt is negative, we still want to return something real so the optimizer 
    # will work. But we need to make sure it returns a *definitely incorrect* but not super 
    # divergent value. This method seems to work, but seems a bit sketchy.
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

    # When optimizing, we are just trying to match the computed eigenfrequencies with the 
    # target eigenfrequencies, so no need to compute the eigenvectors
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

    sortedresults =
        sort([(bcleanedvalues[i], bvectors[:, i]) for i in 1:length(bcleanedvalues)])
    if axis == ẑ
        # axial modes are sorted from lowest to highest eigenfrequency
        return sortedresults
    else
        # this ordering is reversed for radial modes
        return reverse(sortedresults)
    end
end

diagonalize_Kij(
    ions::Vector{<:Ion},
    axis::NamedTuple{(:x, :y, :z)},
    kz::Real,
    P::Real,
    kr::Real;
    optimize=false
) = diagonalize_Kij(
    [mass(ion) for ion in ions],
    [ion.speciesproperties.charge for ion in ions],
    axis,
    kz,
    P,
    kr,
    optimize=optimize
)

#=
This performs the same function as diagonalize_Kij, but for a homogeneous chain where
all of the ions have the same mass. In that case, one can write the Jacobian Kij in
terms of the characteristic_frequencies, which, in this case are the uniqueCOM frequencies.
Then we can just diagonalize Kij directly instead of having to optimize.
=#
function diagonalize_Kij_for_homogeneous_chain(
    N::Int,
    com::NamedTuple{(:x, :y, :z)},
    axis::NamedTuple{(:x, :y, :z)}
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
    M::Vector{<:Real},
    Q::Vector{<:Int},
    target_eigenfrequencies::NamedTuple{(:x, :y, :z)},
    kz::Real,
    P::Real,
    kr::Real,
    just_y=false
)
    @assert length(M) == length(Q) "length(M) ≠ length(Q)"
    N = length(M)
    axes = just_y ? [ŷ] : [x̂, ŷ]
    axes_dict = Dict(x̂ => :x, ŷ => :y)
    computed_eigenfrequencies = []
    for i in eachindex(axes)
        axis = axes[i]
        computed_eigenfrequency = diagonalize_Kij(M, Q, axis, kz, P, kr, optimize=true)
        push!(computed_eigenfrequencies, computed_eigenfrequency)
    end
    difference = 0
    for i in eachindex(computed_eigenfrequencies)
        difference +=
            abs(computed_eigenfrequencies[i] - target_eigenfrequencies[axes_dict[axes[i]]])
    end
    return difference
end

sequentially_diagonalize(
    ions::Vector{<:Ion},
    target_eigenfrequencies::NamedTuple{(:x, :y, :z)},
    kz::Real,
    P::Real,
    kr::Real,
    just_y=false
) = sequentially_diagonalize(
    [mass(ion) for ion in ions],
    [ion.speciesproperties.charge for ion in ions],
    target_eigenfrequencies,
    kz,
    P,
    kr,
    just_y=just_y
)

function searchfor_trappingpotential_parameters(
    M::Vector{<:Real},
    Q::Vector{<:Int},
    target_eigenfrequencies::NamedTuple{(:x, :y, :z)}
)
    @assert length(M) == length(Q) "length(M) ≠ length(Q)"
    maxM = maximum(M)
    value = diagonalize_Kij(M, Q, ẑ, maxM, 0, 0; optimize=true)
    kz = (target_eigenfrequencies.z / value)^2 * maxM
    lowerbounds = [-1, -1]
    upperbounds = [Inf, Inf]
    difference = 0
    if target_eigenfrequencies.x != target_eigenfrequencies.y
        initial = [10kz, kz / 2]
        res = optimize(
            x ->
                sequentially_diagonalize(M, Q, target_eigenfrequencies, kz, x[1], x[2]),
            lowerbounds,
            upperbounds,
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
            lowerbounds[1:1],
            upperbounds[1:1],
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
        "Solver failed to find a normal mode structure compatible with \
        `characteristicfrequencies`."
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
    ions::Vector{<:Ion},
    target_eigenfrequencies::NamedTuple{(:x, :y, :z)}
) = searchfor_trappingpotential_parameters(
    [mass(ion) for ion in ions],
    [ion.speciesproperties.charge for ion in ions],
    target_eigenfrequencies
)


#############################################################################################
# Object fields
#############################################################################################

"""
    full_normal_mode_description(chain::LinearChain)

For each axis, this contains an array of tuples where the first element is a vibrational
frequency [Hz] and the second element is a vector describing the participation of each ion
at that vibrational frequency (i.e. the normal mode eigenvector corresponding to that 
eigenfrequency).
"""
full_normal_mode_description(lc::LinearChain) =
    full_normal_mode_description(ions(lc), comfrequencies(lc))

"""
    full_normal_mode_description(ions, comfreqs::NamedTuple{(:x, :y, :z)})

Same thing but with an iterable list of `ions` and `NamedTuple` of COM frequencies
explicitly given.
"""
function full_normal_mode_description(ions, comfreqs::NamedTuple{(:x, :y, :z)})
    try
        ions = collect(ions)
    catch
        AssertionError("`ions` must be an iterable collection of `<:Ion`.") |> throw
    end
    @assert !isnan(comfreqs.x) "This function doesn't work for user defined mode structure."
    trappingparams = searchfor_trappingpotential_parameters(ions, comfreqs)
    normalmodes = (
        x=diagonalize_Kij(ions, x̂, trappingparams...),
        y=diagonalize_Kij(ions, ŷ, trappingparams...),
        z=diagonalize_Kij(ions, ẑ, trappingparams...)
    )
    return normalmodes
end

"""
    function full_normal_mode_description(
        M::Vector{<:Real}, Q::Vector{<:Int}, comfreqs::NamedTuple{(:x, :y, :z)}
    )
Same thing but explicitly provide the masses `M` and charges `Q` of the ions.
"""
function full_normal_mode_description(
    M::Vector{<:Real},
    Q::Vector{<:Int},
    comfreqs::NamedTuple{(:x, :y, :z)}
)
    @assert length(M) == length(Q) "length(M) ≠ length(Q)"
    trappingparams = searchfor_trappingpotential_parameters(M, Q, comfreqs)
    normalmodes = (
        x=diagonalize_Kij(M, Q, x̂, trappingparams...),
        y=diagonalize_Kij(M, Q, ŷ, trappingparams...),
        z=diagonalize_Kij(M, Q, ẑ, trappingparams...)
    )
    return normalmodes
end

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
function _construct_vibrational_modes(x, num_ions)
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
            el = x[xyz[i]]
            if !isempty(el)
                # allows slices like (y=[1, 2:5], )
                elscreened = []
                for k in eachindex(el)
                    if typeof(el[k]) == Colon
                        @assert length(el) == 1 "Improper indexing."
                        elscreened = 1:num_ions
                    else
                        push!(elscreened, el[k])
                    end
                end
                flattenedx = reduce(
                    vcat,
                    [typeof(i) <: UnitRange ? collect(i) : i for i in elscreened]
                )
            else
                flattenedx = []
            end
            push!(values, flattenedx)
        end
    end
    return (; zip(xyz, values)...)
end

Base.length(lc::LinearChain) = length(lc.ions)

"""
    ionpositions(chain::LinearChain)
Returns the positions of the ions in `chain` in meters.
"""
ionpositions(chain) = chain.ionpositions


#############################################################################################
# Visualization
#############################################################################################

"""
    visualize(lc::LinearChain, axis::String, modes::Vector{<:Int}, format="bars")

Visualize the normal mode structure of a `lc` as a Plots.jl plot. `axis` can either be one
of "x", "y", "z" or a NamedTuple{(:x,:y,:z)}. `modes` is a vector of indices 
selecting which modes should be included in the plot. Slicing is supported. `format` can be 
either "bars" or "circles."

Note that the indexing refers to the full normal mode description and not the subset 
`selectedmodes(lc)`.
"""
function visualize(lc::LinearChain, axis, modes; format="bars", legend=true)
    ion_list = ions(lc)
    axes_dict = Dict("x" => :x, "y" => :y, "z" => :z, x̂ => :x, ŷ => :y, ẑ => :z)
    direction = axes_dict[axis]
    radialplane = direction in [:x, :y] ? true : false
    fnm = full_normal_mode_description(ion_list, comfrequencies(lc))[direction]
    num_ions = length(ion_list)
    unique_ion_labels = unique([ion.speciesproperties.shortname for ion in ion_list])
    color_list = []
    for ion in ion_list
        push!(color_list, findfirst(ion.speciesproperties.shortname .== unique_ion_labels))
    end

    indextypes = map(typeof, modes)
    if Colon in indextypes
        indices = [1:num_ions...]
    elseif UnitRange{Int} in indextypes
        indices = reduce(vcat, [typeof(i) <: UnitRange ? collect(i) : i for i in modes])
    else
        indices = modes
    end
    modes = fnm[indices]
    subplots = []
    legend_plot = plot()
    for label in unique_ion_labels
        plot!(
            [],
            [],
            color_palette=palette(:tab10),
            label=" " * label,
            lw=10,
            legend=:top,
            legendcolumns=1,
            grid=false,
            showaxis=false,
            legendfontsize=15,
            foreground_color_legend=nothing
        )
    end
    for mode in modes
        frequency = round(mode[1] / 1e6, digits=3)
        v = mode[2]
        if format == "bars"
            push!(
                subplots,
                bar(
                    v,
                    xaxis=false,
                    title="$frequency MHz",
                    legend=false,
                    color=color_list,
                    color_palette=palette(:tab10),
                    label=false,
                    yticks=false,
                    showaxis=false
                )
            )
            hline!([0], label=false, lc=:black)
        elseif format == "circles"
            v = map(x -> x == 0 ? NaN : x, v)
            offset = 0.6
            v = map(x -> x + sign(x) * offset, v)
            xpos = linear_equilibrium_positions(num_ions)
            GR.setarrowsize(0.75)
            if !radialplane
                push!(
                    subplots,
                    quiver(
                        xpos,
                        zeros(num_ions),
                        quiver=(v / 3, zeros(num_ions)),
                        lw=2,
                        color=:black
                    )
                )
            else
                push!(
                    subplots,
                    quiver(
                        xpos,
                        zeros(num_ions),
                        quiver=(zeros(num_ions), v),
                        lw=2,
                        color=:black
                    )
                )
            end
            scatter!(
                xpos,
                zeros(length(xpos)),
                markersize=15,
                bg=:white,
                showaxis=false,
                grid=false,
                markerstrokewidth=0,
                ylim=radialplane ? [minimum([v..., 0] .- 0.5), maximum([v..., 0] .+ 0.5)] :
                     [-0.5, 0.5],
                xlim=[xpos[1] - 0.5, xpos[end] + 0.5],
                size=radialplane ? (500, 500) : (1200, 300),
                title="$frequency MHz",
                legend=false,
                color_palette=palette(:tab10),
                markercolor=color_list
            )
        end
    end
    nrows = 2
    divbythree = length(modes) ÷ nrows
    modthree = length(modes) % nrows
    if divbythree == 0
        l = (modthree + 1, 1)
    elseif modthree == 0
        if legend
            l = @layout [grid(1, 1){0.1w} grid(nrows, divbythree)]
        else
            l = @layout [grid(nrows, divbythree)]
        end
    else
        if legend
            l = @layout [grid(1, 1){0.1w} grid(nrows, divbythree) grid(modthree, 1)]
        else
            l = @layout [grid(nrows, divbythree) grid(modthree, 1)]
        end
    end
    if legend
        p = plot(legend_plot, subplots..., layout=l)
    else
        p = plot(subplots..., layout=l)
    end
    display(p)
end

"""
    visualize(vm::VibrationalMode; format="bars")

Visualize the normal mode structure of a `vm` as a Plots.jl plot. `format` can be either
`bar` or `ions`. `format` can be either "bars" or "circles."
"""
function visualize(vm::VibrationalMode; format="bars")
    axes_dict = Dict(x̂ => :x, ŷ => :y, ẑ => :z)
    direction = axes_dict[vm.axis]
    radialplane = direction in [:x, :y] ? true : false
    frequency = round(vm.ν / 1e6, digits=3)
    v = vm.modestructure
    if format == "bars"
        bar(
            v,
            xticks=(1:4, ones(length(v))),
            xaxis=false,
            ylabel="Relative participation\nof ion in mode",
            title="$frequency MHz",
            legend=false
        )
        hline!([0], label=false, lc=:black)
    elseif format == "circles"
        v = map(x -> x == 0 ? NaN : x, v)
        N = length(v)
        xpos = linear_equilibrium_positions(N)
        GR.setarrowsize(0.75)
        scatter(
            xpos,
            zeros(length(xpos)),
            markersize=15,
            bg=:white,
            showaxis=false,
            grid=false,
            markerstrokecolor=:white,
            markercolor=:red3,
            ylim=radialplane ? [minimum([v..., 0] .- 0.5), maximum([v..., 0] .+ 0.5)] :
                 [-0.5, 0.5],
            xlim=[xpos[1] - 0.5, xpos[end] + 0.5],
            size=(1200, 300),
            title="$frequency MHz",
            legend=false,
        )
        offset_vector = copy(v)
        for i in eachindex(offset_vector)
            if offset_vector[i] > 0
                offset_vector[i] += 0.1
            else
                offset_vector[i] -= 0.1
            end
        end
        if !radialplane
            quiver!(
                xpos,
                zeros(N),
                quiver=(offset_vector / 2, zeros(N)),
                lw=2,
                color=:black
            )
        else
            quiver!(xpos, zeros(N), quiver=(zeros(N), offset_vector), lw=3, color=:black)
        end
    end
    display(current())
end

"""
    visualize(
        vm::Tuple{Float64, Vector{Float64}}, axis::NamedTuple{(:x, :y, :z)}; format="bars"
    )

Same thing but input a normal mode description as a tuple with first element the 
eigenfrequency and second the eigenvector.
"""
function visualize(
    vm::Tuple{Float64, Vector{Float64}},
    axis::NamedTuple{(:x, :y, :z)};
    format="bars"
)
    VM = VibrationalMode(vm..., axis=axis, N=1)
    visualize(VM; format=format)
end
