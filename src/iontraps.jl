using LinearAlgebra: eigen
using NLsolve: nlsolve
using .PhysicalConstants: e, ϵ₀

export IonTrap,
    ions,
    linear_equilibrium_positions,
    Anm,
    LinearChain,
    characteristic_length_scale,
    comfrequencies,
    selectedmodes,
    full_normal_mode_description,
    modes,
    xmodes,
    ymodes,
    zmodes,
    modecutoff!

"""
    IonTrap
Physical configuration of ions. Stores a collection of ions and information about the
interactions of their center of mass motion.
"""
abstract type IonTrap end

# required functions
ions(I::IonTrap) = I.ions

#############################################################################################
# LinearChain - a linear Coulomb crystal
#############################################################################################

"""
    linear_equilibrium_positions(N::Int)
Returns the scaled equilibrium positions of `N` ions in a harmonic potential, assuming that
all ions have the same mass.
[ref](https://doi.org/10.1007/s003400050373)
"""
function linear_equilibrium_positions(N::Int)
    function f!(F, x, N)
        for i in 1:N
            F[i] = (
                x[i] - sum([1 / (x[i] - x[j])^2 for j in 1:(i - 1)]) + sum([1 / (x[i] - x[j])^2 for j in (i + 1):N])
            )
        end
    end

    function j!(J, x, N)
        for i in 1:N, j in 1:N
            if i ≡ j
                J[i, j] = (
                    1 + 2 * sum([1 / (x[i] - x[j])^3 for j in 1:(i - 1)]) - 2 * sum([1 / (x[i] - x[j])^3 for j in (i + 1):N])
                )
            else
                J[i, j] = (
                    -2 * sum([1 / (x[i] - x[j])^3 for j in 1:(i - 1)]) + 2 * sum([1 / (x[i] - x[j])^3 for j in (i + 1):N])
                )
            end
        end
    end
    # see eq.8 in the ref to see where (2.018/N^0.559) comes from
    if isodd(N)
        initial_x = [(2.018 / N^0.559) * i for i in (-N ÷ 2):(N ÷ 2)]
    else
        initial_x =
            [(2.018 / N^0.559) * i for i in filter(x -> x ≠ 0, collect((-N ÷ 2):(N ÷ 2)))]
    end
    sol = nlsolve((F, x) -> f!(F, x, N), (J, x) -> j!(J, x, N), initial_x, method = :newton)
    return sol.zero
end

"""
    characteristic_length_scale(M::Real, ν::Real)
Returns the characteristic length scale for a linear chain of identical ions of mass `M`
and with axial trap frequency 2π × ν.
"""
characteristic_length_scale(M::Real, ν::Real) = (e^2 / (4π * ϵ₀ * M * (2π * ν)^2))^(1 / 3)

"""
    Anm(N::Real, com::NamedTuple{(:x,:y,:z)}, axis::NamedTuple{(:x,:y,:z)})
Computes the normal modes and corresponding trap frequencies along a particular `axis` for a
collection of `N` ions in a linear Coloumb crystal and returns an array of tuples with first
element the frequency of the normal mode and 2nd element the corresponding eigenvector.

`com` should be a `NamedTuple` of COM frequences for the different axes:
`(x<:Real, y<:Real, z<:Real)`, where the ``z``-axis is taken to be parallel to the axis of
the crystal.
"""
function Anm(N::Int, com::NamedTuple{(:x, :y, :z)}, axis::NamedTuple{(:x, :y, :z)})
    @assert axis in [x̂, ŷ, ẑ] "axis can be x̂, ŷ, ẑ"
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
        _sparsify!(view(a2, :, i), 1e-5)
    end
    if axis == ẑ
        return [(a1[i], a2[:, i]) for i in 1:length(a1)]
    else
        a = reverse([(a1[i], a2[:, i]) for i in 1:length(a1)])
        return a
    end
end

_sparsify!(x, eps) = @. x[abs(x) < eps] = 0

"""
    LinearChain(;
            ions::Vector{Ion},
            comfrequencies::NamedTuple{(:x,:y,:z)},
            selectedmodes::NamedTuple{(:x,:y,:z),Tuple{Vararg{Vector{VibrationalMode},3}}}
        )

Contains all of the information necessary to describe a collection of ions trapped in a 3D
harmonic potential and forming a linear coulomb crystal.

**user-defined fields**
* `ions::Vector{Ion}`: a list of ions that compose the linear Coulomb crystal
* `comfrequencies::NamedTuple{(:x,:y,:z),Tuple{Vararg{Vector{VibrationalMode},3}}}`:
    Describes the COM frequencies `(x=ν_x, y=ν_y, z=ν_z)`. The ``z``-axis is taken to be
    parallel to the crystal's symmetry axis.
* `selectedmodes::NamedTuple{(:x,:y,:z)}`:  e.g. `selectedmodes=(x=[1], y=[2], z=[1,2])`.
    Specifies the axis and a list of integers which correspond to the ``i^{th}`` farthest
    mode away from the COM for that axis. For example, `selectedmodes=(z=[2])` would
    specify the axial stretch mode. These are the modes that will be modeled in the chain.
    Note: `selectedmodes=(x=[],y=[],z=[1])`, `selectedmodes=(y=[],z=[1])`
    and `selectedmodes=(;z=[1])` are all acceptable and equivalent.
"""
struct LinearChain <: IonTrap  # Note: this is not a mutable struct
    ions::Vector{<:Ion}
    comfrequencies::NamedTuple{(:x, :y, :z)}
    selectedmodes::NamedTuple{(:x, :y, :z), Tuple{Vararg{Vector{VibrationalMode}, 3}}}
    function LinearChain(; ions, comfrequencies::NamedTuple, selectedmodes::NamedTuple, N=10)
        selectedmodes = _construct_vibrational_modes(selectedmodes)
        warn = nothing
        for i in 1:length(ions), j in (i + 1):length(ions)
            @assert typeof(ions[i]) == typeof(ions[j]) "multispecies chains not yet supported; all ions in chain must be of same species"
            if ions[j] ≡ ions[i]
                ions[j] = copy(ions[i])
                if isnothing(warn)
                    warn = "Some ions point to the same thing. Making copies."
                    @warn warn
                end
            end
        end
        nions = length(ions)
        A = (
            x = Anm(nions, comfrequencies, x̂),
            y = Anm(nions, comfrequencies, ŷ),
            z = Anm(nions, comfrequencies, ẑ)
        )
        vm = (
            x = Vector{VibrationalMode}(undef, 0),
            y = Vector{VibrationalMode}(undef, 0),
            z = Vector{VibrationalMode}(undef, 0)
        )
        r = [x̂, ŷ, ẑ]
        for (i, modes) in enumerate(selectedmodes), mode in modes
            push!(vm[i], VibrationalMode(A[i][mode]..., axis = r[i], N=N))
        end
        l = linear_equilibrium_positions(length(ions))
        l0 = characteristic_length_scale(mass(ions[1]), comfrequencies.z) # Needs to be changed when allowing for multi-species chains. Current workaround takes the mass of only the first ion to define the characteristic length scale.
        for (i, ion) in enumerate(ions)
            Core.setproperty!(ion, :ionnumber, i)
            Core.setproperty!(ion, :ionposition, l[i] * l0)
        end
        return new(ions, comfrequencies, vm)
    end
end

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
    full_normal_mode_description(chain::LinearChain)::NamedTuple{(:x,:y,:z)}
For each axis, this contains an array of tuples where the first element is a vibrational
frequency [Hz] and the second element is a vector describing the corresponding normalized
normal mode structure.
"""
function full_normal_mode_description(chain::LinearChain)
    nions = length(ions(chain))
    com_freqs = comfrequencies(chain)
    A = (x = Anm(nions, com_freqs, x̂),
        y = Anm(nions, com_freqs, ŷ),
        z = Anm(nions, com_freqs, ẑ)
    )
    return A
end

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
Returns an array of all of the selected `VibrationalModes` in the x-direction in the `LinearChain`.
"""
xmodes(lc::LinearChain) = selectedmodes(lc).x
"""
    ymodes(lc::LinearChain)
Returns an array of all of the selected `VibrationalModes` in the y-direction in the `LinearChain`.
"""
ymodes(lc::LinearChain) = selectedmodes(lc).y
"""
    zmodes(lc::LinearChain)
Returns an array of all of the selected `VibrationalModes` in the z-direction in the `LinearChain`.
"""
zmodes(lc::LinearChain) = selectedmodes(lc).z

"""
    modecutoff!(lc::LinearChain, N::Int)
Sets the upper bound of the Hilbert space of all `VibrationalMode`s in `lc` to be the Fock state `N`.
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
    @assert isnothing(findfirst(x -> x ∉ xyz, k)) "keys of `selectedmodes` must be `:x`, `:y` or `:z`"
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
