# IonTrap


Once you've defined a number of ions`<:Ion`, the next step is to put them together into an `IonTrap`. 

**DocString**: [`IonTrap`](@ref)

The following subtypes of IonTraps are currently implemented:

+ [`LinearChain`](@ref): A Linear Coulomb crystal, which may be composed of ions that have different mass and or charge.

And the following configurations are under consideration for implementation:

+ *CoulombCrystal*: A two-dimensional ion crystal.
+ *Ring*: A one-dimensional chainof ions in a ring configuration.