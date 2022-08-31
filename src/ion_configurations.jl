using LinearAlgebra: eigen
using NLsolve: nlsolve
using Optim
using .PhysicalConstants: e, ϵ₀

export IonConfiguration,
    ions,
    LinearChain,
    get_vibrational_modes,
    axial_trap_frequency,
    radial_trap_frequency


"""
    IonConfiguration

Physical configuration of ions. Stores a collection of ions and information about the 
interactions of their center of mass motion.
"""
abstract type IonConfiguration end

# required functions
ions(I::IonConfiguration)::Vector{Ion} = I.ions
                                                  

#############################################################################################
# LinearChain - a linear Coulomb crystal
#############################################################################################

#=
    linear_equilibrium_positions(N::Int, M::Vector{<:Real}; withJacobian=false)

Returns the scaled equilibrium positions of `N` ions with masses `M` in a harmonic potential.
To return actual spacings in meters, multiply results by:
```math 
l = e² / 4πε₀Mν²
```
where ``M`` equals the maximum of `M` and ν is the frequency defining the harmonic potential.
[ref](https://doi.org/10.1007/s003400050373).
If `withJacobian`` is true, the analytic Jacobian is used in the solver (but this doesn't 
generally seem to be necessary).
=#
function linear_equilibrium_positions(N::Int; withJacobian=false)
    function f!(F, x, N)
        for i in 1:N
            F[i] = (
                  x[i] 
                - sum([1 / (x[i] - x[j])^2 for j in 1:(i - 1)]) 
                + sum([1 / (x[i] - x[j])^2 for j in (i + 1):N])
            )
        end
    end

    function j!(J, x, N)
        for i in 1:N, j in 1:N
            if i ≡ j
                J[i, j] = (
                      1 
                    + 2 * (sum([1 / (x[i] - x[j])^3 for j in 1:(i - 1)]) 
                    - 2 * sum([1 / (x[i] - x[j])^3 for j in (i + 1):N]))
                )
            else
                J[i, j] = (
                     -2 * (sum([1 / (x[i] - x[j])^3 for j in 1:(i - 1)]) 
                    - 2 * sum([1 / (x[i] - x[j])^3 for j in (i + 1):N]))
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
    if withJacobian
        sol = nlsolve(
            (F, x) -> f!(F, x, N), (J, x) -> j!(J, x, N), initial_x, method = :newton
        )
    else
        sol = nlsolve((F, x) -> f!(F, x, N), initial_x, method = :newton)
    end
    @assert sol.f_converged "Failed to find equilibrium positions for linear chain."
    return sol.zero
end

#=
Computes the eigenvalues/vectors for the Hessian, undoes the mass weighting and returns in 
sorted order (Ascending in terms of eigenvalue if axial, otherwise descending. If the ions 
are all the same mass, this ordering gives special significance to the center-of-mass modes 
-- which are the lowest frequency in the axial direction and highest in the radial).

+ M: A vector of the ion masses.

+ com: A named tuple with entries which correspond to the com frequencies of the ion chain 
*only* when the ions are all the same mass. Otherwise this corresponds to the lowest 
eigenfrequency in the axial direction and the highest eigenfrequency in the radial directions.

+ axis: The axis along which the Hessian is computed (q = x, y or z).

+ ω: A vector of the ion vibrational frequencies along the relevant trap axis. This is the 
natural frequency of vibration that the ion would experience along the relevant axis if it 
was in the trap alone. If this is equal to nothing, then these values are set internally 
according to psuedopotential_constant and DC_imbalance.

+ psuedopotential_constant =  

+ DC_imbalance = 

+ k_axial: The (mass-independent) spring constant for the axial modes.
=#
function linear_chain_normal_modes(
    M::Vector{<:Real}, 
    com::Union{NamedTuple{(:x, :y, :z)}, Nothing}, 
    axis::NamedTuple{(:x, :y, :z)};
    ω::Union{Vector, Nothing}=nothing,
    psuedopotential_constant=0, 
    DC_imbalance=1,
    k_axial=1
)
    # Constuct Hessian H (https://doi.org/10.1007/s100530170275)
    N = length(M)
    l = linear_equilibrium_positions(N)
    if isnothing(com)
        k_axial = maximum(M) * (2π * 1e6)^2
    end
    a, β = axis == ẑ ? (2, 0) : (-1, psuedopotential_constant * maximum(M))
    H = Array{Real}(undef, N, N)
    for n in 1:N, j in 1:N
        if n ≡ j
            if isnothing(ω)
                H[n, j] = (
                    (-1 / 2)^(axis != ẑ) * DC_imbalance + β / M[n] + 
                    a * sum([1 / abs(l[j] - l[p])^3 for p in 1:N if p != j])
                )
            else
                H[n, j] = (
                    M[n] * ω[n]^2 / k_axial + 
                    a * sum([1 / abs(l[j] - l[p])^3 for p in 1:N if p != j])
                )
            end
        else
            H[n, j] = -a / abs(l[j] - l[n])^3
        end
        # Mass-weight the matrix elements. This means the eigenvectors 
        # of this Hessian will need to be 'unweighted' later
        H[n, j] /= √(M[n] * M[j])
    end
    h = eigen(H)
    h1 = []
    for (i, value) in enumerate(h.values)
        if axis == ẑ
            @assert i > 0 "Outside of linear chain stability regime (negative eigenvalues)"
            push!(h1, sqrt(value))
        else
            push!(h1, value > 0 ? sqrt(value) : -sqrt(-value))
        end
    end
    if axis == ẑ && !isnothing(com)
        # the axial spring constants are mass-independent
        k_axial = (2π * com.z / minimum(h1))^2 # k_axial: axial spring constant
    end
    h1 *= √k_axial / 2π
    h2 = h.vectors
    for i in 1:size(h2, 2)
        for j in 1:length(h2[:, i])
            h2[j, i] /= √M[j]
        end
        h2[:, i] .= normalize(h2[:, i])
    end
    if axis == ẑ
        return k_axial, sort([(h1[i], h2[:, i]) for i in 1:length(h1)])
    else
        return reverse(sort([(h1[i], h2[:, i]) for i in 1:length(h1)]))
    end
end

function minimize(f, axis, com)
    for i in 1:1e4
        res = optimize(
            x -> abs(f(x)[1][1] - com[axis]), 
            [i * randn()], 
            Optim.Options(g_abstol=1e-12, g_reltol=1e-12)
        )
        (res.iterations == 0) && continue
        !res.g_converged && continue
        a = f(res.minimizer)
        for ai in a
            if ai[1] < 0
                error("Outside of linear chain stability regime (negative eigenvalues)")
            end
        end
        return res.minimizer[1]
    end
    return false
end

function minimize_psuedopotential_constant(M, com, k_axial)
    if all(y -> y == M[1], M)
        pc = nothing
    else
        pc = minimize(
            x -> linear_chain_normal_modes(
                    M, com, x̂, k_axial=k_axial, 
                    psuedopotential_constant=first(x), DC_imbalance=1
                ),
            :x, com
        )
    end
    if isnothing(pc)
        dci = -2 * (com.x / com.z)^2
        pc = 0
    else
        dci = 1
    end
    a = linear_chain_normal_modes(
            M, com, x̂, k_axial=k_axial, psuedopotential_constant=pc, DC_imbalance=dci
        )
    return pc, dci, a
end

function minimize_DC_imbalance(M, com, k_axial, pc)
    if all(y -> y == M[1], M)
        dci = nothing
    else
        dci = minimize(
            x -> linear_chain_normal_modes(
                    M, com, ŷ, k_axial=k_axial, 
                    psuedopotential_constant=pc, DC_imbalance=first(x)
                ),
            :y, com
        )
    end
    if isnothing(dci)
        dci = -2 * (com.y / com.z)^2
        pc = 0
    end
    a = linear_chain_normal_modes(
                M, com, ŷ, k_axial=k_axial, psuedopotential_constant=pc, DC_imbalance=dci
            )
    return pc, dci, a
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
        Describes the COM frequencies `(x=ω_x, y=ω_y, z=ω_z)`. The ``z``-axis is taken to be 
        parallel to the crystal's symmetry axis and we assume (but don't directly enforce) 
        that ``ω\\_z/ω\\_{x,y} > 0.73N^{0.86}`` [ref](https://doi.org/10.1007/s003400050225)
        When the ions are different species (so different masses) there is no longer generally 
        a unique com mode. In this case we take com to correspond to the lowest (highest) 
        eigenfrequency of the axial (radial) modes, which corresponds to the homogoneous case.
* `vibrational_modes::NamedTuple{(:x,:y,:z)}`:  eg. `(x=[1], y=[2], z=[1,2])`. 
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
    ion_frequencies::NamedTuple{(:x, :y, :z), <:Tuple} 
    vibrational_modes::NamedTuple{(:x, :y, :z), Tuple{Vararg{Vector{VibrationalMode}, 3}}}
    full_normal_mode_description::NamedTuple{(:x, :y, :z)}
    function LinearChain(; 
            ions, com_frequencies=nothing, ion_frequencies=nothing, 
            vibrational_modes::NamedTuple
        )
        vibrational_modes = _construct_vibrational_modes(vibrational_modes)
        warn = nothing
        for i in 1:length(ions), j in (i + 1):length(ions)
            if ions[j] ≡ ions[i]
                ions[j] = copy(ions[i])
                if isnothing(warn)
                    # only output warning once
                    warn = "Some ions point to the same thing. Making copies."
                    @warn warn
                end
            end
        end
        M = mass.(ions)
        if isnothing(ion_frequencies)
            k_axial, z = linear_chain_normal_modes(M, com_frequencies, ẑ)
            pcx, dcix, x = minimize_psuedopotential_constant(M, com_frequencies, k_axial)  
            pcy, dciy, y = minimize_DC_imbalance(M, com_frequencies, k_axial, pcx)
            xf, yf, zf = [], [], []
            for ionmass in M
                ωx = sqrt(k_axial * (-dcix / 2 + pcx * maximum(M) / ionmass) / ionmass)
                ωy = sqrt(k_axial * (-dciy / 2 + pcy * maximum(M) / ionmass) / ionmass)
                ωz = sqrt(k_axial / ionmass)
                push!(xf, ωx)
                push!(yf, ωy)
                push!(zf, ωz)
            end
            ion_frequencies = (x=xf, y=yf, z=zf)
        else
            k_axial, z = linear_chain_normal_modes(M, nothing, ẑ, ω=ion_frequencies.z)
            x = linear_chain_normal_modes(M, nothing, x̂, ω=ion_frequencies.x, k_axial=k_axial)
            y = linear_chain_normal_modes(M, nothing, ŷ, ω=ion_frequencies.y, k_axial=k_axial)
            com_frequencies = (x=maximum(first.(x)), y=maximum(first.(y)), z=minimum(first.(z)))
        end
        A = (x=x, y=y, z=z)
        vm = (
            x = Vector{VibrationalMode}(undef, 0),
            y = Vector{VibrationalMode}(undef, 0),
            z = Vector{VibrationalMode}(undef, 0)
        )
        r = [x̂, ŷ, ẑ]
        for (i, modes) in enumerate(vibrational_modes), mode in modes
            push!(vm[i], VibrationalMode(A[i][mode]..., axis = r[i]))
        end
        l = linear_equilibrium_positions(length(ions))
        l0 = characteristic_length_scale(k_axial)
        for (i, ion) in enumerate(ions)
            Core.setproperty!(ion, :ionnumber, i)
            Core.setproperty!(ion, :position, l[i] * l0)
        end
        return new(ions, com_frequencies, ion_frequencies, vm, A)
    end
end

"""
    get_vibrational_modes(lc::LinearChain)

Returns an array of all of the selected `VibrationalModes` in the `LinearChain`.
The order is `[lc.x..., lc.y..., lc.z...]`.
"""
get_vibrational_modes(lc::LinearChain) = collect(Iterators.flatten(lc.vibrational_modes))

function Base.print(lc::LinearChain)
    println("$(length(lc.ions)) ions")
    println("com frequencies: $(lc.com_frequencies)")
    println("selected vibrational_modes: $(lc.vibrational_modes)")
end

# suppress long output
Base.show(io::IO, lc::LinearChain) = print(io, "LinearChain($(length(lc.ions)) ions)")

# Takes e.g. (y=[1]) to (x=[], y=[1], z=[])
function _construct_vibrational_modes(x)
    k = collect(keys(x))
    xyz = [:x, :y, :z]
    @assert isnothing(findfirst(x -> x ∉ xyz, k)) ("keys of `vibrational_modes` must be " *
    "`:x`, `:y` or `:z`")
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


#############################################################################################
# General functions
#############################################################################################

"""
    axial_trap_frequency(VDC::Real, M::Real, z₀:Real, ξ::Real)

Computes the ion trap frequency along the axial direction according to:

``ω_z = \\sqrt{2 e V_{DC} ξ / M z₀²}``

where `e` is the charge of the electron, `VDC` is the voltage on the endcap electrodes, 
`M` is the mass of the ion, `z₀` is the distance from the endcap 
electrodes to the trapping point and `ξ` is a dimensionless constant between 0 and 1 that 
accounts for deviations from an ideal trap. [ref](https://doi.org/10.1063/1.367318)
"""
axial_trap_frequency(VDC, M, z₀, ξ) = √((2 * e * VDC * ξ) / (M * z₀^2))

"""
    radial_trap_frequency(
        VDC::Real, VRF::Real, ΩRF::Real, M::Real, z₀::Real, r₀::Real, ξ::Real, ψ::Real
    )

Computes the ion trap frequency along one of the radial directions according to:

``ω_{x,y} = \\sqrt{ω_{DC}^2/2 + (ψeV_{RF} / 2Mr₀²Ω_{RF})²}``

``ω_{DC} = `` `axial_trap_frequency(VDC, M, z₀, ξ)`

where `e` is the charge of the electron, `VDC` is the voltage on the endcap electrodes, 
`VRF` is the amplitude of the voltage on the RF electrodes, `ΩRF` is (2π ×) the RF 
frequency, `z₀` is the distance from the endcap electrodes to the trapping point, `r₀` is the 
distance from the RF electrodes to the trapping point `ξ` is a dimensionless constant between 
0 and 1 that accounts for deviations of the DC fields from an ideal trap and `ψ` accounts for 
similar deviations for the RF fields. [ref](https://doi.org/10.1063/1.367318)
"""
function radial_trap_frequency(VDC, VRF, ΩRF, M, z₀, r₀, ξ, ψ) 
    ωDC = axial_trap_frequency(VDC, M, z₀, ξ)
    return √((-ωDC^2 / 2 + (ψ * e * VRF / (2 * M * r₀^2 * ΩRF))^2))
end

#=
    characteristic_length_scale(M::Real, ν::Real)

Returns the characteristic length scale for a linear chain of identical ions of mass `M`
and with axial trap frequency ``2π × ν``.
=#
characteristic_length_scale(M::Real, ν::Real) = (e^2 / (4π * ϵ₀ * M * (2π * ν)^2))^(1 / 3)
#=
    characteristic_length_scale(k::Real) 
    
``k = Mν²`` 
=#
characteristic_length_scale(k::Real) = (e^2 / (4π * ϵ₀ * k))^(1 / 3)
