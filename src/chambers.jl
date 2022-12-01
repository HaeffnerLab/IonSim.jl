using QuantumOptics: tensor, CompositeBasis
using .PhysicalConstants: ħ, c

export Chamber,
    basis,
    iontrap,
    lasers,
    ions,
    modes,
    xmodes,
    ymodes,
    zmodes,
    modecutoff!,
    groundstate,
    Efield_from_pi_time,
    Efield_from_pi_time!,
    transition_frequency,
    wavelength!,
    wavelength_from_transition!,
    set_gradient!,
    Efield_from_rabi_frequency,
    Efield_from_rabi_frequency!,
    global_beam!,
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
        `[T.iontrap.vibrational_modes.x, T.iontrap.vibrational_modes.y,
        T.iontrap.vibrational_modes.z]`.
    E.g. for:

    ```
    chain = LinearChain(ions=[C1, C2], com_frequencies=(x=2e6,y=2e6,z=1e6),
    selected_modes=(x=[1, 2], y=[], z=[1]))
    ```

    The ordering of the basis would be

    `C1.basis ⊗ C2.basis ⊗ chain.vibrational_modes.x[1].basis
    ⊗ chain.vibrational_modes.x[2].basis ⊗ chain.vibrational_modes.z[1].basis`

    Otherwise, the ordering is according to the form of the initial state used in the solver.
"""
mutable struct Chamber
    iontrap::LinearChain
    B::Real
    Bhat::NamedTuple{(:x, :y, :z)}
    ∇B::Real
    δB::Function
    lasers::Array{<:Laser}
    basis::CompositeBasis
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
        basis = tensor(
            iontrap.ions...,
            iontrap.vibrational_modes.x...,
            iontrap.vibrational_modes.y...,
            iontrap.vibrational_modes.z...,
        )
        return new(iontrap, B, Bhat, ∇B, δBt, lasers, basis, _cnst_δB)
    end
end

function Base.setproperty!(T::Chamber, s::Symbol, v::Tv) where {Tv}
    if s == :δB
        if Tv <: Number
            _cnst_δB = true
            vt(t) = v
        else
            _cnst_δB = false
            vt = v
        end
        Core.setproperty!(T, s, vt)
        Core.setproperty!(T, :_cnst_δB, _cnst_δB)
    elseif s == :basis || s == :_cnst_δB
        return
    elseif s == :iontrap
        Core.setproperty!(T, s, v)
        basis = tensor(
            v.ions...,
            v.vibrational_modes.x...,
            v.vibrational_modes.y...,
            v.vibrational_modes.z...,
        )
        Core.setproperty!(T, :basis, basis)
    else
        Core.setproperty!(T, s, v)
    end
end

Base.show(io::IO, T::Chamber) = print(io, "Chamber")  # suppress long output


#############################################################################################
# Object fields
#############################################################################################
iontrap(T::Chamber) = T.iontrap
lasers(T::Chamber) = T.lasers


#############################################################################################
# general functions
#############################################################################################

"""
    ionintrap(trap::Chamber, ion::Ion)
Returns a boolean that indicates whether `ion` is actually in `trap`. Useful for checking if an error needs to be thrown.
"""
function ionintrap(trap::Chamber, ion::Ion)
    return ion in ions(trap.iontrap)
end

"""
    basis(T::Chamber)
Returns the composite basis describing the Hilbert space for `T`.
"""
function basis(T::Chamber)::CompositeBasis # Isn't this already constructed as T.basis?
    return tensor(
        T.iontrap.ions...,
        T.iontrap.vibrational_modes.x...,
        T.iontrap.vibrational_modes.y...,
        T.iontrap.vibrational_modes.z...,
    )
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
If obj is a `LinearChain` or `Chamber`, returns the full ground state of the motional degrees of freedom as a tensor product.
"""
groundstate(mode::VibrationalMode) = mode[0]
groundstate(modes::Vector{VibrationalMode}) = tensor([mode[0] for mode in modes]...)
groundstate(lc::LinearChain) = groundstate(modes(lc))
groundstate(T::Chamber) = groundstate(modes(T))


"""
    global_beam!(T::Chamber, laser::Laser)
Set `laser` to shine with full intensity on all ions in `Chamber`.
"""
function global_beam!(T::Chamber, laser::Laser)
    for n in eachindex(T.iontrap.ions)
        push!(laser.pointing, (n, 1.0))
    end
end

"""
    Efield_from_pi_time(
        pi_time::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion,
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
Compute the E-field needed to get a certain `pi_time` with a certain `laser`-`ion`
`transition`.

Alternatively, one may use
```
Efield_from_pi_time(
            pi_time::Real, T::Chamber, laser_index::Int, ion_index::Int,
            transition::Union{Tuple{String,String},Vector{<:String}}
        )
```
which is the same as
`Efield_from_pi_time(pi_time, T.Bhat, T.lasers[laser_index], T.iontrap.ions[ion_index],
transition)`
"""
function Efield_from_pi_time(
    pi_time::Real,
    Bhat::NamedTuple{(:x, :y, :z)},
    laser::Laser,
    ion::Ion,
    transition::Tuple
)
    p = laser.pointing
    s_indx = findall(x -> x[1] == ionnumber(ion), p)
    @assert length(s_indx) > 0 "This laser doesn't shine on this ion"
    s = p[s_indx[1]][2]
    Ω = s * matrix_element(ion, transition, 1, laser.k, laser.ϵ, Bhat)
    if Ω < 1e-15
        # even when coupling strength is zero, numerical error causes it to be finite
        # (on order 1e-16), this is a band-aid to prevent users from unknowingly setting
        # the E-field to something absurd (like 1e20 V/m)
        return Inf
    end
    return 1 / (2Ω * pi_time)
end

function Efield_from_pi_time(
    pi_time::Real,
    T::Chamber,
    laser_index::Int,
    ion_index::Int,
    transition::Tuple
)
    return Efield_from_pi_time(
        pi_time,
        T.Bhat,
        T.lasers[laser_index],
        T.iontrap.ions[ion_index],
        transition
    )
end

"""
    Efield_from_pi_time!(
            pi_time::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion,
            transition::Union{Tuple{String,String},Vector{<:String}}
        )
Same as [`Efield_from_pi_time`](@ref), but updates `laser[:E]` in-place.
"""
function Efield_from_pi_time!(
    pi_time::Real,
    Bhat::NamedTuple{(:x, :y, :z)},
    laser::Laser,
    ion::Ion,
    transition::Tuple
)
    Efield::Float64 = Efield_from_pi_time(pi_time, Bhat, laser, ion, transition)
    return laser.E = t -> Efield
end

function Efield_from_pi_time!(
    pi_time::Real,
    T::Chamber,
    laser_index::Int,
    ion_index::Int,
    transition::Tuple
)
    Efield::Float64 = Efield_from_pi_time(pi_time, T, laser_index, ion_index, transition)
    return T.lasers[laser_index].E = t -> Efield
end

"""
    Efield_from_rabi_frequency(
        Ω::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion,
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
Compute the Efield needed to get a certain rabi frequency `Ω` with a certain `laser`-`ion`
`transition`.

Alternatively, one may use
```
Efield_from_rabi_frequency(
            Ω::Real, T::Chamber, laser_index::Int, ion_index::Int,
            transition::Union{Tuple{String,String},Vector{<:String}}
        )
```
which is the same as
`Efield_from_rabi_frequency(pi_time, T.Bhat, T.lasers[laser_index],
T.iontrap.ions[ion_index], transition)`
"""
function Efield_from_rabi_frequency(
    Ω::Real,
    Bhat::NamedTuple{(:x, :y, :z)},
    laser::Laser,
    ion::Ion,
    transition::Tuple
)
    return Efield_from_pi_time(1 / 2Ω, Bhat, laser, ion, transition)
end

function Efield_from_rabi_frequency(
    Ω::Real,
    T::Chamber,
    laser_index::Int,
    ion_index::Int,
    transition::Tuple
)
    return Efield_from_rabi_frequency(
        Ω,
        T.Bhat,
        T.lasers[laser_index],
        T.iontrap.ions[ion_index],
        transition
    )
end

"""
    Efield_from_rabi_frequency!(
        Ω::Real, Bhat::NamedTuple{(:x,:y,:z)}, laser::Laser, ion::Ion,
        transition::Union{Tuple{String,String},Vector{<:String}}
    )
Same as [`Efield_from_rabi_frequency`](@ref), but updates `laser[:E]` in-place.
"""
function Efield_from_rabi_frequency!(
    Ω::Real,
    Bhat::NamedTuple{(:x, :y, :z)},
    laser::Laser,
    ion::Ion,
    transition::Tuple
)
    Efield::Float64 = Efield_from_rabi_frequency(Ω, Bhat, laser, ion, transition)
    return laser.E = t -> Efield
end

function Efield_from_rabi_frequency!(
    Ω::Real,
    T::Chamber,
    laser_index::Int,
    ion_index::Int,
    transition::Tuple
)
    Efield::Float64 = Efield_from_rabi_frequency(Ω, T, laser_index, ion_index, transition)
    return T.lasers[laser_index].E = t -> Efield
end

"""
    Bfield(T::Chamber, ion::Ion)
Retuns the value of the magnetic field in `T` at the location of `ion`, including both the trap's overall B-field and its B-field gradient.
"""
function Bfield(T::Chamber, ion::Ion)
    @assert ionintrap(T, ion) "trap does not contain ion"
    return T.B + T.∇B * ionposition(ion)
end

"""
    transitionfrequency(ion::Ion, transition::Tuple, T::Chamber; ignore_manualshift=false)
Retuns The frequency of the transition `transition` including the Zeeman shift experienced by `ion` given its position in `T`.

One may alternatively replace `ion` with `ion_index::Int`, which instead specifies the index of the intended ion within `T`.
"""
transitionfrequency(ion::Ion, transition::Tuple, T::Chamber; ignore_manualshift = false) =
    transitionfrequency(
        ion,
        transition;
        B = Bfield(T, ion),
        ignore_manualshift = ignore_manualshift
    )
transitionfrequency(ion_index::Int, transition::Tuple, T::Chamber; ignore_manualshift = false) =
    transitionfrequency(
        T.iontrap.ions[ion_index],
        transition,
        T;
        ignore_manualshift = ignore_manualshift
    )

"""
    transitionwavelength(ion::Ion, transition::Tuple, T::Chamber; ignore_manualshift=false)
Retuns The wavelength of the transition `transition` including the Zeeman shift experienced by `ion` given its position in `T`.

One may alternatively replace `ion` with `ion_index`::Int, which instead specifies the index of the intended ion within `T`.
"""
transitionwavelength(ion::Ion, transition::Tuple, T::Chamber; ignore_manualshift = false) =
    transitionwavelength(
        ion,
        transition;
        B = Bfield(T, ion),
        ignore_manualshift = ignore_manualshift
    )
transitionwavelength(
    ion_index::Int,
    transition::Tuple,
    T::Chamber;
    ignore_manualshift = false
) = transitionwavelength(
    T.iontrap.ions[ion_index],
    transition,
    T;
    ignore_manualshift = ignore_manualshift
)

"""
    wavelength!(laser::Laser, wavelength::Real)
Sets the wavelength of `laser` to `wavelength`.
"""
function wavelength!(laser::Laser, wavelength::Real)
    laser.λ = wavelength
end

"""
    wavelength_from_transition!(laser::Laser, ion::Ion, transition::Tuple, Bfield::Real)
    wavelength_from_transition!(laser::Laser, ion::Ion, transition::Tuple, T::Chamber)
Sets the wavelength of `laser` to the transition wavelength of `transition` in the ion `ion`,
at a magnetic field value given by `Bfield` if the final argument is a `Real`, or at the magnetic
field seen by `ion` if the final argument is the `Chamber` that contains it.
"""
function wavelength_from_transition!(laser::Laser, ion::Ion, transition::Tuple, Bfield::Real)
    wavelength = transitionwavelength(ion, transition, B=Bfield)
    laser.λ = wavelength
    return wavelength
end
function wavelength_from_transition!(laser::Laser, ion::Ion, transition::Tuple, T::Chamber)
    wavelength = transitionwavelength(ion, transition, T)
    laser.λ = wavelength
    return wavelength
end

"""
    matrix_element(I::Ion, transition::Tuple, T::Chamber, laser::Laser, time::Real)
Calls `matrix_element(I::Ion, transition::Tuple, Efield::Real, khat::NamedTuple, ϵhat::NamedTuple, Bhat::NamedTuple=(;z=1))`
with `Efield`, `khat`, and `ϵhat` evaluated for `laser` at time `time`, and `Bhat` evaluated for `T`.

One may alternatively replace `ion` with `ion_index`::Int, which instead specifies the index of the intended ion within `T`.
"""
matrix_element(I::Ion, transition::Tuple, T::Chamber, laser::Laser, time::Real) =
    matrix_element(I, transition, laser.E(time), laser.k, laser.ϵ, T.Bhat)
matrix_element(ion_index::Int, transition::Tuple, T::Chamber, laser::Laser, time::Real) =
    matrix_element(T.iontrap.ions[ion_index], transition, T, laser, time)

"""
    zeeman_shift(I::Ion, sublevel, T::Chamber)
Calls `zeeman_shift(I::Ion, sublevel, B::Real)` with `B` evaluated for ion `I` in `T`.

One may alternatively replace `ion` with `ion_index`::Int, which instead specifies the index of the intended ion within `T`.
"""
zeeman_shift(I::Ion, sublevel::Union{Tuple{String, Real}, String}, T::Chamber) =
    zeeman_shift(I, sublevel, Bfield(T, I))
zeeman_shift(ion_index::Int, sublevel::Union{Tuple{String, Real}, String}, T::Chamber) =
    zeeman_shift(T.iontrap.ions[ion_index], sublevel, T)

"""
    set_gradient!(
            T::Chamber, ion_indxs::Tuple{Int,Int}, transition::Tuple, f::Real
        )
Sets the Bfield gradient in place to achieve a detuning `f` between the `transition` of two
ions, which are assumed to be of the same species. `ion_indxs` refer to the
ordering of the ions in the chain.
"""
function set_gradient!(T::Chamber, ion_indxs::Tuple{Int, Int}, transition::Tuple, f::Real)
    ionA = T.iontrap.ions[ion_indxs[1]]
    ionB = T.iontrap.ions[ion_indxs[2]]
    separation = abs(ionposition(ionA) - ionposition(ionB))

    (SL1, SL2) = transition
    L1 = sublevel2level(ionA, SL1)
    L2 = sublevel2level(ionA, SL2)
    g1 = landegf(ionA, L1)
    g2 = landegf(ionA, L2)
    m1 = quantumnumbers(ionA, SL1).m
    m2 = quantumnumbers(ionA, SL2).m
    # Calculate Zeeman shifts with a unit B-field using a method of zeeman_shift that ensures a nonlinear term is not used
    E1 = zeeman_shift(1.0, g1, m1)
    E2 = zeeman_shift(1.0, g2, m2)
    return T.∇B = f / (abs(E2 - E1) * separation)
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
    lambdicke(V::VibrationalMode, L::Laser, I::Ion)
The Lamb-Dicke parameter: 
``|k|cos(\\theta)\\sqrt{\\frac{\\hbar}{2m\\nu}}`` 
for a given vibrational mode, ion and laser.
"""
function lambdicke(V::VibrationalMode, L::Laser, I::Ion; scaled = false)
    @fastmath begin
        k = 2π / L.λ
        scaled ? ν = 1 : ν = V.ν
        x0 = √(ħ / (2 * mass(I) * 2π * ν))
        cosθ = ndot(L.k, V.axis)
        k * x0 * cosθ * V.mode_structure[ionnumber(I)]
    end
end
