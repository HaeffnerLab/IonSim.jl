using LinearAlgebra: eigen
using NLsolve: nlsolve
using .PhysicalConstants: e, ϵ₀

export IonConfiguration,
    ions,
    linear_equilibrium_positions,
    Anm,
    LinearChain,
    characteristic_length_scale,
    get_vibrational_modes

"""
    IonConfiguration
Physical configuration of ions. Stores a collection of ions and  information about the
interactions of their center of mass motion.
"""
abstract type IonConfiguration end

# required functions
ions(I::IonConfiguration)::Vector{Ion} = I.ions

#############################################################################################
# LinearChain - a linear Coulomb crystal
#############################################################################################

"""
    linear_equilibrium_positions(N::Int)
Returns the scaled equilibrium positions of `N` ions in a harmonic potential.
If ions don't all have the same mass,
specify a list of masses (or relative masses) for each ion.

If the mass list is shorter than the number of ions, it is assumed to repeat;
i.e. (M_a, M_b) corresponds to (M_a, M_b, M_a, M_b, ...).
[ref](https://doi.org/10.1007/s003400050373)
"""
function linear_equilibrium_positions(N::Int, masslist::NTuple{num, Real} = (1.0,)) where {num}
    relativemasslist = masslist./masslist[1]
    masslistlength = length(relativemasslist)
    function f!(F, x, N, masses, masseslength)
        for i in 1:N
            F[i] = (
                masses[(i-1) % masseslength + 1] * x[i] - sum([1 / (x[i] - x[j])^2 for j in 1:(i - 1)]) + sum([1 / (x[i] - x[j])^2 for j in (i + 1):N])
            )
        end
    end

    function j!(J, x, N, masses, masseslength)
        for i in 1:N, j in 1:N
            if i ≡ j
                J[i, j] = (
                    masses[(i-1) % masseslength + 1] + 2 * sum([1 / (x[i] - x[j])^3 for j in 1:(i - 1)]) - 2 * sum([1 / (x[i] - x[j])^3 for j in (i + 1):N])
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
    sol = nlsolve((F, x) -> f!(F, x, N, relativemasslist, masslistlength),
                  (J, x) -> j!(J, x, N, relativemasslist, masslistlength),
                   initial_x, method = :newton)
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
    # what should change here?
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
            ions::Vector{Ion}, com_frequencies::NamedTuple{(:x,:y,:z)},
            vibrational_modes::NamedTuple{(:x,:y,:z),Tuple{Vararg{Vector{VibrationalMode},3}}}
        )

Contains all of the information necessary to describe a collection of ions trapped in a 3D
harmonic potential and forming a linear coulomb crystal.

**user-defined fields**
* `ions::Vector{Ion}`: a list of ions that compose the linear Coulomb crystal
* `com_frequencies::NamedTuple{(:x,:y,:z),Tuple{Vararg{Vector{VibrationalMode},3}}}`:
        Describes the COM frequencies `(x=ν_x, y=ν_y, z=ν_z)`. The ``z``-axis is taken to be
        parallel to the crystal's symmetry axis.
* `vibrational_modes::NamedTuple{(:x,:y,:z)}`:  e.g. `vibrational_modes=(x=[1], y=[2], z=[1,2])`.
    Specifies the axis and a list of integers which correspond to the ``i^{th}`` farthest
    mode away from the COM for that axis. For example, `vibrational_modes=(z=[2])` would
    specify the axial stretch mode. These are the modes that will be modeled in the chain.
    Note: `vibrational_modes=(x=[],y=[],z=[1])`, `vibrational_modes=(y=[],z=[1])`
    and `vibrational_modes=(;z=[1])` are all acceptable and equivalent.
**derived fields**
* `full_normal_mode_description::NamedTuple{(:x,:y,:z)}`: For each axis, this contains an
    array of tuples where the first element is a vibrational frequency [Hz] and the second
    element is a vector describing the corresponding normalized normal mode structure.
"""
struct LinearChain <: IonConfiguration  # Note: this is not a mutable struct
    ions::Vector{<:Ion}
    com_frequencies::NamedTuple{(:x, :y, :z)}
    vibrational_modes::NamedTuple{(:x, :y, :z), Tuple{Vararg{Vector{VibrationalMode}, 3}}}
    full_normal_mode_description::NamedTuple{(:x, :y, :z)}
    function LinearChain(; ions, com_frequencies, vibrational_modes::NamedTuple)
        vibrational_modes = _construct_vibrational_modes(vibrational_modes)
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
        N = length(ions)
        A = (
            x = Anm(N, com_frequencies, x̂),
            y = Anm(N, com_frequencies, ŷ),
            z = Anm(N, com_frequencies, ẑ)
        )
        vm = (
            x = Vector{VibrationalMode}(undef, 0),
            y = Vector{VibrationalMode}(undef, 0),
            z = Vector{VibrationalMode}(undef, 0)
        )
        r = [x̂, ŷ, ẑ]
        for (i, modes) in enumerate(vibrational_modes), mode in modes
            push!(vm[i], VibrationalMode(A[i][mode]..., axis = r[i]))
        end
        l = linear_equilibrium_positions(length(ions), Tuple(mass.(ions)))
        # Use the mass of only the first ion to define the characteristic length scale.
        l0 = characteristic_length_scale(mass(ions[1]), com_frequencies.z)
        for (i, ion) in enumerate(ions)
            Core.setproperty!(ion, :ionnumber, i)
            Core.setproperty!(ion, :position, l[i] * l0)
        end
        return new(ions, com_frequencies, vm, A)
    end
end

"""
    get_vibrational_modes(lc::LinearChain)
Returns an array of all of the selected `VibrationalModes` in the `LinearChain`.
The order is `[lc.x..., lc.y..., lc.z...]`.
"""
function get_vibrational_modes(lc::LinearChain)
    return collect(Iterators.flatten(lc.vibrational_modes))
end

function Base.print(lc::LinearChain)
    println("$(length(lc.ions)) ions")
    println("com frequencies: $(lc.com_frequencies)")
    return println("selected vibrational_modes: $(lc.vibrational_modes)")
end

function Base.show(io::IO, lc::LinearChain)  # suppress long output
    return print(io, "LinearChain($(length(lc.ions)) ions)")
end

# Takes e.g. (y=[1]) to (x=[], y=[1], z=[])
function _construct_vibrational_modes(x)
    k = collect(keys(x))
    xyz = [:x, :y, :z]
    @assert isnothing(findfirst(x -> x ∉ xyz, k)) "keys of `vibrational_modes` must be `:x`, `:y` or `:z`"
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
