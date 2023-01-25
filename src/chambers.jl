using QuantumOptics: tensor, CompositeBasis
using .PhysicalConstants: ħ, c

export Chamber,
    iontrap,
    bfield,
    bfield_unitvector,
    bgradient,
    bfield_fluctuation,
    lasers,
    basis,
    iontrap!,
    bfield!,
    bfield_unitvector!,
    bgradient!,
    bfield_fluctuation!,
    lasers!,
    ions,
    modes,
    xmodes,
    ymodes,
    zmodes,
    modecutoff!,
    groundstate,
    intensity_from_pitime,
    intensity_from_pitime!,
    intensity_from_rabifrequency,
    intensity_from_rabifrequency!,
    transition_frequency,
    wavelength_from_transition!,
    globalbeam!,
    lambdicke

#############################################################################################
# an ion trap with a linear ion chain configuration
#############################################################################################

"""
    Chamber(;
            iontrap::LinearChain, B::Real=0, Bhat::NamedTuple{(:x,:y,:z)=ẑ, ∇B::Real=0,
             δB::Union{Real,Function}=0, lasers::Vector{Laser}
    )

Information necessary to describe the Hamiltonian for a collection of ions in a linear chain
interacting with laser light.
**user-defined fields**
* `iontrap<:LinearChain`
* `B`: A real value describing the mean magnitude of the B-field [Tesla].
* `Bhat::NamedTuple{(:x,:y,:z)}`: Describes the direction of the B-field (defaults to ẑ).
* `∇B`: Magnitude of the B-field gradient. We assume that the gradient always points along the
        z-direction. [Tesla / meter]
* `δB::Function`: Time-dependence of the B-field [Tesla]
* `lasers::Array{<:Laser}`: For each laser in the array, the pointing field should contain
        an array of `Tuple{Int,Real}`. The first element specifies the index of an ion
        in the `ions` field that the laser interacts with. The second element specifies a
        scaling factor for the strength of that interaction (to be used, e.g., for
        modeling cross-talk).
**derived fields**
* `_cnst_δB::Bool`: A Boolean flag signifying whether or not `δB` is a constant function.
* `basis<:CompositeBasis`: The basis for describing the combined system, ions + vibrational
        modes. If constructing the Hamiltonian explictly (with [`hamiltonian`](@ref)), then
        the ordering of the basis is set, by convention, as
        ``ion₁ ⊗ ion₂ ⊗ ... ⊗ ion_N ⊗ mode₁ ⊗ mode₂ ⊗ ... ⊗ mode_N``, where the ion bases are
        ordered according to the order in `T.iontrap.ions` and the vibrational modes
        are ordered according to the order in
        `[T.iontrap.selectedmodes.x, T.iontrap.selectedmodes.y,
        T.iontrap.selectedmodes.z]`.
    E.g. for:

    ```
    chain = LinearChain(ions=[C1, C2], comfrequencies=(x=2e6,y=2e6,z=1e6),
    selectedmodes=(x=[1, 2], y=[], z=[1]))
    ```

    The ordering of the basis would be

    `C1.basis ⊗ C2.basis ⊗ chain.selectedmodes.x[1].basis
    ⊗ chain.selectedmodes.x[2].basis ⊗ chain.selectedmodes.z[1].basis`

    Otherwise, the ordering is according to the form of the initial state used in the solver.
"""
mutable struct Chamber
    iontrap::LinearChain
    B::Real
    Bhat::NamedTuple{(:x, :y, :z)}
    ∇B::Real
    δB::Function
    lasers::Array{<:Laser}
    _cnst_δB::Bool
    function Chamber(;
        iontrap::LinearChain,
        B = 0,
        Bhat = ẑ,
        ∇B = 0,
        δB::TδB = 0,
        lasers = Laser[]
    ) where {TδB}
        warn = nothing
        for i in 1:length(lasers)
            if length(lasers[i].pointing) == 0
                for n in eachindex(iontrap.ions)
                    push!(lasers[i].pointing, (n, 1.0))
                end
            end
            for j in (i + 1):length(lasers)
                if lasers[j] ≡ lasers[i]
                    lasers[j] = copy(lasers[i])
                    if isnothing(warn)
                        warn = "Some lasers point to the same thing. Making copies."
                        @warn warn
                    end
                end
            end
        end
        @assert isapprox(norm(Bhat), 1, rtol = 1e-6) "!(|$Bhat| = 1)"
        for (li, l) in enumerate(lasers), p in l.pointing
            @assert p[1] <= length(iontrap.ions) (
                """lasers[$li] points at iontrap.ions[$(p[1])], but there are only
                 $(length(iontrap.ions)) ions."""
            )
        end
        if TδB <: Number
            _cnst_δB = true
            δBt(t) = δB
        else
            _cnst_δB = false
            δBt = δB
        end
        return new(iontrap, B, Bhat, ∇B, δBt, lasers, _cnst_δB)
    end
end

Base.show(io::IO, T::Chamber) = print(io, "Chamber")  # suppress long output


#############################################################################################
# Object fields
#############################################################################################

iontrap(chamber::Chamber) = chamber.iontrap
bfield(chamber::Chamber) = chamber.B
bfield_unitvector(chamber::Chamber) = chamber.Bhat
bgradient(chamber::Chamber) = chamber.∇B
bfield_fluctuation(chamber::Chamber) = chamber.δB
lasers(chamber::Chamber) = chamber.lasers


#############################################################################################
# Setters
#############################################################################################

function iontrap!(chamber::Chamber, iontrap::IonTrap)
    chamber.iontrap = iontrap
end

function bfield!(chamber::Chamber, B::Real)
    chamber.B = B
end

function bfield_unitvector!(chamber::Chamber, Bhat::NamedTuple{(:x, :y, :z)})
    rtol = 1e-6
    @assert isapprox(norm(Bhat), 1, rtol = rtol) "!(|̂B| = 1)"
    chamber.Bhat = Bhat
end

function bgradient!(chamber::Chamber, ∇B::Real)
    chamber.∇B = ∇B
end

function bfield_fluctuation!(chamber::Chamber, δB::Function)
    chamber.δB = δB
    chamber._cnst_δB = false
end
function bfield_fluctuation!(chamber::Chamber, δB::Real)
    chamber.δB = (t -> δB)
    chamber._cnst_δB = true
end

function lasers!(chamber::Chamber, lasers::Vector{Laser})
    chamber.lasers = lasers
end


#############################################################################################
# general functions
#############################################################################################


"""	
    basis(chamber::Chamber)	
Returns the composite basis describing the Hilbert space for `chamber`.
This is the same as basis(iontrap(chain)).
"""	
function basis(T::Chamber)
    return tensor(	
        T.iontrap.ions...,	
        T.iontrap.selectedmodes.x...,	
        T.iontrap.selectedmodes.y...,	
        T.iontrap.selectedmodes.z...,	
    )	
end


"""
    ionintrap(trap::Chamber, ion::Ion)
Returns a boolean that indicates whether `ion` is actually in `trap`. Useful for checking if an error needs to be thrown.
"""
function ionintrap(trap::Chamber, ion::Ion)
    return ion in ions(trap.iontrap)
end

""""
    ions(T::Chamber)
Returns a list of the ions in the `Chamber`.
"""
ions(T::Chamber) = ions(iontrap(T))


"""
    modes(T::Chamber)
Returns modes(iontrap(T))
"""
modes(T::Chamber) = modes(iontrap(T))
"""
    xmodes(T::Chamber)
Returns an array of all of the selected `VibrationalModes` in the x-direction in the `Chamber`'s `IonConfiguration`.
"""
xmodes(T::Chamber) = xmodes(iontrap(T))
"""
    ymodes(T::Chamber)
Returns an array of all of the selected `VibrationalModes` in the y-direction in the `Chamber`'s `IonConfiguration`.
"""
ymodes(T::Chamber) = ymodes(iontrap(T))
"""
    zmodes(T::Chamber)
Returns an array of all of the selected `VibrationalModes` in the z-direction in the `Chamber`'s `IonConfiguration`.
"""
zmodes(T::Chamber) = zmodes(iontrap(T))

"""
    modecutoff!(T::Chamber, N::Int)
Sets the upper bound of the Hilbert space of all `VibrationalMode`s in the `IonTrap` of `T` to be the Fock state `N`.
"""
function modecutoff!(T::Chamber, N::Int)
    modecutoff!(iontrap(T), N)
end

"""
    groundstate(obj)
If obj is a `VibrationalMode`, returns the N=0 ket of that mode.
If obj is a Vector of `VibrationalMode`, returns a tensor product `mode1[0] ⊗ mode2[0] ⊗ ...` in the same order given.
If obj is a `LinearChain`, returns the full ground state of the motional degrees of freedom as a tensor product.
"""
groundstate(mode::VibrationalMode) = mode[0]
groundstate(modes::Vector{VibrationalMode}) = tensor([mode[0] for mode in modes]...)
groundstate(lc::LinearChain) = groundstate(modes(lc))


"""
    globalbeam!(laser, chamber::Chamber)
Set `laser` to shine with full intensity on all ions in `Chamber`.
`laser` may be either a Laser or an Int indicating the desired laser's index within `chamber`.
"""
function globalbeam!(laser::Laser, chamber::Chamber)
    p = [(n, 1.0) for n in eachindex(ions(chamber))]
    pointing!(laser, p)
end
function globalbeam!(laserindex::Int, chamber::Chamber)
    laser = lasers(chamber)[laserindex]
    globalbeam!(laser, chamber)
end


"""
    intensity_from_pitime(
        laser::Laser, pi_time::Real, ion::Ion, transition::Tuple,
        Bhat::NamedTuple{(:x,:y,:z)}
    )
Compute the intensity needed to get a certain `pi_time` with a certain resonant `laser`-`ion`
`transition`, in the presence of a magnetic field pointing in the direction `Bhat`.
"""
function intensity_from_pitime(
    laser::Laser,
    pi_time::Real,
    ion::Ion,
    transition::Tuple,
    Bhat::NamedTuple{(:x, :y, :z)}
)
    p = pointing(laser)
    s_indx = findall(x -> x[1] == ionnumber(ion), p)
    @assert length(s_indx) > 0 "This laser doesn't shine on this ion"
    s = p[s_indx[1]][2]
    Ω = s * matrixelement(ion, transition, 1.0, polarization(laser), wavevector(laser), Bhat)
    if Ω < 3e-14    # After change from Efield to intensity: This inequality changed so that it serves the same function but now for intensity = 1.0 rather than efield = 1.0 (specificially, increased by a factor of √(2/(cϵ₀)) ∼ 30 in SI units)
        # even when coupling strength is zero, numerical error causes it to be finite
        # (on order 1e-14), this is a band-aid to prevent users from unknowingly setting
        # the intensity to something absurd (like 1e20 V/m)
        return Inf
    end
    return (1 / (2Ω * pi_time))^2
end

"""
    intensity_from_pitime(
        laser, pi_time::Real, ion, transition::Tuple, chamber::Chamber
        )
Compute the intensity needed to get a certain `pi_time` with a certain resonant `laser`-`ion`
`transition` within `chamber`, which defines the magnetic field direction.
`laser` may be either a Laser or an Int indicating the desired laser's index within `chamber`.
`ion` may be either an Ion or an Int indicating the desired ion's index within `chamber`.
"""
function intensity_from_pitime(
    laser::Laser,
    pi_time::Real,
    ion::Ion,
    transition::Tuple,
    chamber::Chamber
)
    return intensity_from_pitime(laser, pi_time, ion, transition, bfield_unitvector(chamber))
end
function intensity_from_pitime(
    laser_index::Int,
    pi_time::Real,
    ion_index::Int,
    transition::Tuple,
    chamber::Chamber
)
    return intensity_from_pitime(
        lasers(chamber)[laser_index],
        pi_time,
        ions(chamber)[ion_index],
        transition,
        bfield_unitvector(chamber))
end


"""
    intensity_from_pitime!(
        laser::Laser, pi_time::Real, ion::Ion, transition::Tuple,
        Bhat::NamedTuple{(:x,:y,:z)}
    )
    intensity_from_pitime!(
        laser, pi_time::Real, ion, transition::Tuple, chamber::Chamber
    )
Same as `intensity_from_pitime`, but updates `laser[:I]` in-place.
"""
function intensity_from_pitime!(
    laser::Laser,
    pi_time::Real,
    ion::Ion,
    transition::Tuple,
    Bhat::NamedTuple{(:x, :y, :z)}
)
    I::Float64 = intensity_from_pitime(laser, pi_time, ion, transition, Bhat)
    intensity!(laser, I)
    return I
end
function intensity_from_pitime!(
    laser::Laser,
    pi_time::Real,
    ion::Ion,
    transition::Tuple,
    chamber::Chamber
)
    I::Float64 = intensity_from_pitime(laser, pi_time, ion, transition, chamber)
    intensity!(laser, I)
    return I
end
function intensity_from_pitime!(
    laser_index::Int,
    pi_time::Real,
    ion_index::Int,
    transition::Tuple,
    chamber::Chamber
)
    I::Float64 = intensity_from_pitime(laser_index, pi_time, ion_index, transition, chamber)
    laser = lasers(chamber)[laser_index]
    intensity!(laser, I)
    return I
end

##############################################################################################

"""
    intensity_from_rabifrequency(
        laser::Laser, rabi_frequency::Real, ion::Ion, transition::Tuple,
        Bhat::NamedTuple{(:x,:y,:z)}
    )
Compute the intensity needed to get a certain `rabi_frequency` with a certain resonant `laser`-`ion`
`transition`, in the presence of a magnetic field pointing in the direction `Bhat`.
"""
function intensity_from_rabifrequency(
    laser::Laser,
    rabi_frequency::Real,
    ion::Ion,
    transition::Tuple,
    Bhat::NamedTuple{(:x, :y, :z)}
)
    return intensity_from_pitime(laser, 1/(2*rabi_frequency), ion, transition, Bhat)
end

"""
intensity_from_rabifrequency(
        laser, rabi_frequency::Real, ion, transition::Tuple, chamber::Chamber
        )
Compute the intensity needed to get a certain `rabi_frequency` with a certain resonant `laser`-`ion`
`transition` within `chamber`, which defines the magnetic field direction.
`laser` may be either a Laser or an Int indicating the desired laser's index within `chamber`.
`ion` may be either an Ion or an Int indicating the desired ion's index within `chamber`.
`laser` and `ion` must either both be indices or both their respective Structs.
"""
function intensity_from_rabifrequency(
    laser::Laser,
    rabi_frequency::Real,
    ion::Ion,
    transition::Tuple,
    chamber::Chamber
)
    return intensity_from_rabifrequency(laser, rabi_frequency, ion, transition, bfield_unitvector(chamber))
end
function intensity_from_rabifrequency(
    laser_index::Int,
    rabi_frequency::Real,
    ion_index::Int,
    transition::Tuple,
    chamber::Chamber
)
    return intensity_from_rabifrequency(
        lasers(chamber)[laser_index],
        rabi_frequency,
        ions(chamber)[ion_index],
        transition,
        bfield_unitvector(chamber))
end


"""
    intensity_from_rabifrequency!(
        laser::Laser, rabi_frequency::Real, ion::Ion, transition::Tuple,
        Bhat::NamedTuple{(:x,:y,:z)}
    )
    intensity_from_rabifrequency!(
        laser, rabi_frequency::Real, ion, transition::Tuple, chamber::Chamber
    )
Same as `intensity_from_rabifrequency!`, but updates `laser[:I]` in-place.
"""
function intensity_from_rabifrequency!(
    laser::Laser,
    rabi_frequency::Real,
    ion::Ion,
    transition::Tuple,
    Bhat::NamedTuple{(:x, :y, :z)}
)
    I::Float64 = intensity_from_rabifrequency(laser, rabi_frequency, ion, transition, Bhat)
    intensity!(laser, I)
    return I
end
function intensity_from_rabifrequency!(
    laser::Laser,
    rabi_frequency::Real,
    ion::Ion,
    transition::Tuple,
    chamber::Chamber
)
    I::Float64 = intensity_from_rabifrequency(laser, rabi_frequency, ion, transition, chamber)
    intensity!(laser, I)
    return I
end
function intensity_from_rabifrequency!(
    laser_index::Int,
    rabi_frequency::Real,
    ion_index::Int,
    transition::Tuple,
    chamber::Chamber
)
    I::Float64 = intensity_from_rabifrequency(laser_index, rabi_frequency, ion_index, transition, chamber)
    laser = lasers(chamber)[laser_index]
    intensity!(laser, I)
    return I
end


"""
    bfield(chamber::Chamber, ion)
Retuns the value of the magnetic field in `T` at the location of `ion`, including both the trap's overall B-field and its B-field gradient.
`ion` may be either an Ion or an Int indicating the desired ion's index within `chamber`.
"""
function bfield(chamber::Chamber, ion::Ion)
    @assert ionintrap(chamber, ion) "trap does not contain ion"
    return bfield(chamber) + bgradient(chamber) * ionposition(ion)
end
function bfield(chamber::Chamber, ion_index::Int)
    return bfield(chamber, ions(chamber)[ion_index])
end

"""
    transitionfrequency(ion, transition::Tuple, chamber::Chamber; ignore_manualshift=false)
Returns The frequency of the transition `transition` including the Zeeman shift experienced by `ion` given its position in `T`.
`ion` may be either an Ion or an Int indicating the desired ion's index within `chamber`.
"""
transitionfrequency(ion::Ion, transition::Tuple, chamber::Chamber; ignore_manualshift = false) =
    transitionfrequency(
        ion,
        transition;
        B = bfield(chamber, ion),
        ignore_manualshift = ignore_manualshift
    )
transitionfrequency(ion_index::Int, transition::Tuple, chamber::Chamber; ignore_manualshift = false) =
    transitionfrequency(
        ions(chamber)[ion_index],
        transition,
        chamber,
        ignore_manualshift = ignore_manualshift
    )

"""
    transitionwavelength(ion, transition::Tuple, chamber::Chamber; ignore_manualshift=false)
Returns The wavelength of the transition `transition` including the Zeeman shift experienced by `ion` given its position in `T`.
`ion` may be either an Ion or an Int indicating the desired ion's index within `chamber`.
"""
transitionwavelength(ion::Ion, transition::Tuple, chamber::Chamber; ignore_manualshift = false) =
    transitionwavelength(
        ion,
        transition;
        B = bfield(chamber, ion),
        ignore_manualshift = ignore_manualshift
    )
transitionwavelength(
    ion_index::Int,
    transition::Tuple,
    chamber::Chamber;
    ignore_manualshift = false
) = transitionwavelength(
    ions(chamber)[ion_index],
    transition,
    chamber;
    ignore_manualshift = ignore_manualshift
)

"""
    wavelength_from_transition!(laser::Laser, ion::Ion, transition::Tuple, B::Real)
Sets the wavelength of `laser` to the transition wavelength of `transition` in the ion `ion`,
at a magnetic field value given by `B`.
"""
function wavelength_from_transition!(laser::Laser, ion::Ion, transition::Tuple, B::Real)
    λ = transitionwavelength(ion, transition, B=B)
    wavelength!(laser, λ)
    return λ
end

"""
    wavelength_from_transition!(laser::Laser, ion, transition::Tuple, chamber::Chamber)
Sets the wavelength of `laser` to the transition wavelength of `transition` in the ion `ion`,
at the magnetic field seen by `ion` in `chamber`.
`ion` may be either an Ion or an Int indicating the desired ion's index within `chamber`.
"""
function wavelength_from_transition!(laser::Laser, ion::Ion, transition::Tuple, chamber::Chamber)
    λ = transitionwavelength(ion, transition, chamber)
    wavelength!(laser, λ)
    return λ
end
function wavelength_from_transition!(laser::Laser, ion_index::Int, transition::Tuple, chamber::Chamber)
    ion = ions(chamber)[ion_index]
    λ = transitionwavelength(ion, transition, chamber)
    wavelength!(laser, λ)
    return λ
end

"""
    matrixelement(ion, transition::Tuple, laser, chamber::Chamber, time::Real)
Calls `matrixelement(ion, transition, I, ϵhat, khat, Bhat)` with `I`, `ϵhat`, and
`khat` evaluated for `laser` at time `time`, and `Bhat` evaluated for `chamber`.

`ion` may be either an Ion or an Int indicating the desired ion's index within `chamber`.
`laser` may be either a Laser or an Int indicating the desired laser's index within `chamber`.
`ion` and `laser` must either both be indices or both their respective Structs.
"""
matrixelement(ion::Ion, transition::Tuple, laser::Laser, chamber::Chamber, time::Real) =
    matrixelement(ion,
        transition,
        intensity(laser)(time),
        polarization(laser),
        wavevector(laser),
        bfield_unitvector(chamber)
    )
matrixelement(ion_index::Int, transition::Tuple, laser_index::Int, T::Chamber, time::Real) =
    matrixelement(ions(T)[ion_index], transition, lasers(T)[laser_index], chamber, time)

"""
    zeemanshift(I, sublevel, T::Chamber)
Calls `zeemanshift(I::Ion, sublevel, B::Real)` with `B` evaluated for ion `I` in `T`.
`I` may be either an Ion or an Int indicating the desired ion's index within `T`.
"""
zeemanshift(I::Ion, sublevel::Union{Tuple{String, Real}, String}, T::Chamber) =
    zeemanshift(I, sublevel, bfield(T, I))
zeemanshift(ion_index::Int, sublevel::Union{Tuple{String, Real}, String}, T::Chamber) =
    zeemanshift(T.iontrap.ions[ion_index], sublevel, T)

"""
    bgradient!(
            T::Chamber, ion_indxs::Tuple{Int,Int}, transition::Tuple, d::Real
        )
Sets the Bfield gradient in place to achieve a detuning `d` between the `transition` of two
ions, which are assumed to be of the same species. `ion_indxs` refer to the
ordering of the ions in the chain.
"""
function bgradient!(T::Chamber, ion_indxs::Tuple{Int, Int}, transition::Tuple, d::Real)
    ionA = ions(T)[ion_indxs[1]]
    ionB = ions(T)[ion_indxs[2]]
    separation = abs(ionposition(ionA) - ionposition(ionB))

    (SL1, SL2) = transition
    L1 = level(ionA, SL1)
    L2 = level(ionA, SL2)
    g1 = landegf(ionA, L1)
    g2 = landegf(ionA, L2)
    m1 = quantumnumbers(ionA, SL1).m
    m2 = quantumnumbers(ionA, SL2).m
    # Calculate Zeeman shifts with a unit B-field using a method of zeemanshift that ensures a nonlinear term is not used
    E1 = zeemanshift(1.0, g1, m1)
    E2 = zeemanshift(1.0, g2, m2)
    bgradient!(T, d / (abs(E2 - E1) * separation))
end

# In QunatumOptics.jl, this method will return true whenever the shapes of b1 and b2 match,
# but we'd like to distinguish, i.e., between Ion ⊗ mode1 ⊗ mode2 and Ion ⊗ mode2 ⊗ mode1
# when mode1.N == mode2.N but mode1.axis ≠ mode2.axis.
function (
    QuantumOptics.:(==)(
        b1::T,
        b2::T
    ) where {T <: CompositeBasis{<:Vector{Int}, <:Tuple{Vararg{<:IonSimBasis}}}}
)
    N = length(b1.bases)
    if N ≠ length(b2.bases)
        return false
    end
    for i in 1:N
        if !(b1.bases[i] == b2.bases[i])
            return false
        end
    end
    return true
end

"""
    lambdicke(V::VibrationalMode, I::Ion, L::Laser)
The Lamb-Dicke parameter: 
``|k|cos(\\theta)\\sqrt{\\frac{\\hbar}{2m\\nu}}`` 
for a given vibrational mode, ion and laser.
"""
function lambdicke(V::VibrationalMode, I::Ion, L::Laser; scaled = false)
    @fastmath begin
        k = 2π / L.λ
        scaled ? ν = 1 : ν = V.ν
        x0 = √(ħ / (2 * mass(I) * 2π * ν))
        cosθ = ndot(wavevector(L), axis(V))
        k * x0 * cosθ * modestructure(V)[ionnumber(I)]
    end
end
