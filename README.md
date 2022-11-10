<p align="center">
  <img src="https://github.com/HaeffnerLab/IonSim.jl/blob/media/logo3_SM.svg?raw=true", width="450px">
</p>

[![test status](https://github.com/HaeffnerLab/IonSim.jl/actions/workflows/test.yml/badge.svg)](https://github.com/HaeffnerLab/IonSim.jl/actions/workflows/test.yml)
[![codecov][codecov-badge]][codecov-url]
[![License: MIT][license-badge]][license-url]
This branch implements the ability for a user to perform an arbitrary rotating frame transformation in the process of generating the hamiltonian. These changes are contained within hamiltonians.jl and do not affect any other files at this time.

This is done by adding an argument to the exported function hamiltonian(), which is called u (for "unitary"). This should be a vector of vectors of tuples. The tuples contain an eigenstate of the bare atomic Hamiltonian at the first index, and an angular frequency (don't forget factors of 2π) at the second index, so that u[n][k][1] gives the kth rotating state of the nth ion, and the angular frequency of the rotation is u[n][k][2].

Inside of hamiltonians.jl, this change is implemented via the replacement of the funcion \_Δmatrix with the similarly-named \_ΔUmatrix. Whereas \_Δmatrix was calculating the detuning of atomic transitions with respect to lasers, \_ΔUmatrix instead returns the quantity u[n][k][2] - u[n][j][2] + νm for the nth ion's transition between eigenstates j and k actuated by the mth laser, and νm is the mth laser's angular frequency. Note that for the interaction picture, u[n][k][2] - u[n][j][2] = νm-2π(Enk - Enj) which is exactly what \_Δmatrix was computing before. Replacing \_Δmatrix with \_ΔUmatrix adequately implements the rotating frame transformation specified by u for all of the off-diagonal elements of the Hamiltonian.

Additionally, a new function is exported called hamiltonian\_interaction\_picture() which simply calls hamiltonian() but does not require a unitary transformation vector as an argument. Instead, it automatically generates the interaction picture transformation and feeds it to hamiltonian().
\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*
A simple tool, built on top of [QuantumOptics.jl](https://qojulia.org/), for simulating the dynamics of a configuration of
trapped ions interacting with laser light.

**IonSim.jl** primarily performs two jobs:
1. Keeps track of the physical parameters necessary for describing the system.
2. Using these parameters, constructs a function that quickly computes the system's Hamiltonian as a function of time.

The functional form of the Hamiltonian can then be used as input to any of the solvers implemented in
[`QuantumOptics.timeevolution`](https://qojulia.org/documentation/timeevolution/timeevolution/). For more information see:

+ Main code: [https://github.com/HaeffnerLab/IonSim.jl/tree/master/src](https://github.com/HaeffnerLab/IonSim.jl/tree/master/src)
+ Documentation: [https://docs.ionsim.org](https://docs.ionsim.org)
+ Examples: [https://examples.ionsim.org](https://examples.ionsim.org)
+ Benchmarks: (coming soon)

If you have an idea for how to improve IonSim, need some help getting things working or have any other IonSim-related questions feel free to open a GitHub issue.

## How to use

1. Install Julia ([instructions here](https://julialang.org/downloads/)).
2. Open up a Terminal session and fire up the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/#The-Julia-REPL-1) with
```bash
$ julia
```
(*If using Windows it's easier to start a REPL with the julia executable available after installation.*)

3. Now run
```julia
julia> using Pkg
julia> Pkg.add("IonSim")
julia> Pkg.add("IJulia")
```
The last line adds Jupyter compatibility, so you'll have the option to start a Julia kernel from Jupyter notebook.

### Updating IonSim

IonSim.jl is a work in progress. You can update your local code to the most recent version
with the following:

```julia
julia> using Pkg
julia> Pkg.update("IonSim")
```

### Extra Linux notes
* Extract the downloaded file
* copy to `/opt` with:

```bash
$ sudo cp -r julia-1.3.1 /opt/
```
(replacing `julia-1.3.1` with the appropriate version)
* then create a symbolic link in `/usr/local/bin` with
```bash
$ sudo ln -s /opt/julia-1.3.1/bin/julia /usr/local/bin/julia
```

## Development

If you want to run IonSim locally:
* Open up the Julia REPL
```bash
$ julia
```
* Open the package manager by pressing `]`
* Run the following:
```julia
pkg> dev IonSim
```
This will clone a version of this repo in `~/.julia/dev/IonSim/`. Then when you make changes to that repo, it will be immediately reflected when using Julia.

You can run tests with
```julia
pkg> test IonSim
```

Don't forget to format the code! `./scripts/format.sh`

To go back to the version in the registry, use
```julia
pkg> free IonSim
```

In order to update the IonSim.jl version that lives in the Julia general registry, change the version number in the Project.toml file, commit the changes and then add a comment to the commit that reads:
```
@JuliaRegistrator register
```

If you have any questions, please make a GitHub issue.

[license-url]: https://github.com/HaeffnerLab/IonSim.jl/blob/master/LICENSE.md
[license-badge]: https://img.shields.io/badge/License-MIT-green.svg

[codecov-url]: https://codecov.io/gh/HaeffnerLab/IonSim.jl
[codecov-badge]: https://codecov.io/gh/HaeffnerLab/IonSim.jl/branch/master/graph/badge.svg

[twitter-url]: https://twitter.com/Berkeley_ions
[twitter-badge]: https://img.shields.io/twitter/follow/Berkeley_ions.svg?style=social&label=@Berkeley_ions

[logo-url]: https://github.com/HaeffnerLab/IonSim.jl/blob/media/smallest_logo.png?raw=true

