# VibrationalMode

**DocsString**: [`VibrationalMode`](@ref)

This object is primarily intended to be constructed during the instantiations of an `IonTrap`.

For example, 

```@example lc2
using IonSim # hide
c = Ca40()
chain = LinearChain(
    ions=[c, c],
    comfrequencies=(x=3e6,y=3e6,z=1e6), 
    selectedmodes=(;x=[1]) 
)

move = xmodes(chain)[1]

typeof(mode) <: VibrationalMode
```

It has the following fields:


* `ν::Real`: frequency [Hz]
* `modestructure::Vector{Real}`: The normalized eigenvector describing the collective motion 
        of the ions belonging to this mode.
* `δν::Union{Function,Real}`: Either a function describing time-dependent fluctuations of `ν`
        or a real number which will be converted to the constant function `t -> δν`.
* `N::Int`: Highest level included in the Hilbert space.
* `axis::NamedTuple{(:x,:y,:z)}`: The axis of symmetry for the vibration. This must lie along
        one of the basis vectors `x̂`, `ŷ` or `ẑ`.
* `shape::Vector{Int}`: Indicates dimension of used Hilbert space (`=[N+1]`).

And can also be utilized to construct quantum states:

```@example lc2
ψ = c ⊗ mode
ψ.bases
```

```@example lc2
ψ.shape
```

```@example lc2
typeof(ψ)
```
