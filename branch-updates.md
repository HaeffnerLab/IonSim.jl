# Ions.jl

+ Reorganized the code in ions.jl.
    - This was done for (my subjective opinion of) readability (with go-ahead from Neil). I tried to stick to the format (from top of file to bottom) 
        1) object fields
        2) object definitions 
        3) helper functions for object definitions 
        4) object modifiers
        5) user-exposed functions that take new objects as arguments
        6) base overrides
    But we should probably split this file up at some point.
+ `full_level_structure` renamed `level_structure` everywhere
+ `full_transitions` renamed `transitions` everywhere
+ Removed `default_sublevels`
+ Removed `Ion.sublevel_aliases`
+ `Ion.sublevels` -> `Ion.selected_sublevels` changing this property now triggers
`_construct_sublevels` to be rerun (though this isn't necessary currently) and updates
`Ion.shape`. Trying to change `Ion.shape` directly triggers an error.
+ New file energylevels.jl that contains a new type `EnergyLevel` and subtypes 
correspondint to different types of atomic states. Currently implemented is 
`LS` for LS-coupled and `J₁K` for J₁K-coupled states. String macros exist to 
build energy levels from standard spectroscopic notation eg.
```julia
julia> ls"4²S_1/2"
|n=4, l=0, s=1/2, j=1/2⟩

julia> jk"¹[5/2]_5/2"
|k=5/2, s=0, j=5/2, f=1⟩
``` 
Or, for hyperfine levels:
```julia
julia> ls"4²S_1/2(F=1)"
|n=4, l=0, s=1/2, j=1/2, f=1⟩

julia> jk"¹[5/2]_5/2(F=1)"
|k=5/2, s=0, j=5/2, f=1⟩
``` 
+ 

+ Remove `ionnumber` and `ionposition` fields from `IonInstance`. Now when an ion is added
to an `IonTrap` an instance of the new type `TrappedIon{S} <: IonInstance{S} where {S<: Any}`
is created that has fields `ionnumber`, which is the same as before, and `iontrap`, which
points to the instance of `IonTrap` that created it. 