# IonSim.jl 

A simple tool, built on top of [QuantumOptics.jl](https://qojulia.org/), for simulating the dynamics of a configuration of 
trapped ions interacting with laser light.

**IonSim.jl** primarily performs two jobs
1. Keeps track of the physical parameters necessary for describing the system.
2. Using these parameters, constructs a function that quickly computes the system's Hamiltonian as a function of time. 

The functional form of the Hamiltonian can then be used as input to any of the solvers implemented in 
[`QuantumOptics.timeevolution`](https://qojulia.org/documentation/timeevolution/timeevolution/). For more information see:

+ The main code: []()
+ Documentation: []()
+ Examples: []()
+ Benchmarks: []()
 
## Installation

### Installing [Julia](https://julialang.org/)
[Platform specific instructions found here](https://julialang.org/downloads/)

### Installing IonSim

Once Julia has been installed, open a terminal session and begin a 
[Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/#The-Julia-REPL-1) with:

```
$ julia
```

And then enter the following commands:

```julia
julia> Using Pkg

julia> Pkg.clone()

julia> Pkg.develop(PackageSpec(path="/path/to/IonSim.jl/"))
```

