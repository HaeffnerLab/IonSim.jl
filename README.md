<p align="center">
  <img src="https://github.com/HaeffnerLab/IonSim.jl/blob/media/logo3_SM.svg?raw=true", width="450px">
</p>

[![Build Status][travis-badge]][travis-url]
[![codecov][codecov-badge]][codecov-url]
[![License: MIT][license-badge]][license-url]
[![Twitter][twitter-badge]][twitter-url]

A simple tool, built on top of [QuantumOptics.jl](https://qojulia.org/), for simulating the dynamics of a configuration of 
trapped ions interacting with laser light.

**IonSim.jl** primarily performs two jobs
1. Keeps track of the physical parameters necessary for describing the system.
2. Using these parameters, constructs a function that quickly computes the system's Hamiltonian as a function of time. 

The functional form of the Hamiltonian can then be used as input to any of the solvers implemented in 
[`QuantumOptics.timeevolution`](https://qojulia.org/documentation/timeevolution/timeevolution/). For more information see:

+ Main code: [https://github.com/HaeffnerLab/IonSim.jl/tree/master/src](https://github.com/HaeffnerLab/IonSim.jl/tree/master/src)
+ Documentation: [https://docs.ionsim.org](https://docs.ionsim.org)
+ Examples: [https://examples.ionsim.org](https://examples.ionsim.org)
+ Benchmarks: 

If you have an idea for how to improve IonSim, need some help getting things working or have any other IonSim-related questions feel free to open up a discussion on our Gitter Channel: [![Gitter](https://badges.gitter.im/IonSim-jl/community.svg)](https://gitter.im/IonSim-jl/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
 
## Installation

### Installing [Julia](https://julialang.org/)
[Platform specific instructions found here](https://julialang.org/downloads/)

If using Linux, once you've extracted the downloaded file copy it to `/opt` with: 

```bash
$ sudo cp -r julia-1.3.1 /opt/
```

(replacing `julia-1.3.1` with the appropriate version) and then create a symbolic link in `/usr/local/bin` with

```bash 
$ sudo ln -s /opt/julia-1.3.1/bin/julia /usr/local/bin/julia
```

### Installing IonSim

Once Julia has been installed, open a terminal session and begin a 
[Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/#The-Julia-REPL-1) with:

```
$ julia
```

(*If using Windows it's easier to start a REPL with the julia executable available after 
installation.*)


And then enter the following commands:

```julia
julia> using Pkg

julia> Pkg.add(PackageSpec(url="https://github.com/HaeffnerLab/IonSim.jl.git"))
```

You'll also need to download the QuantumOptics package:

```julia
julia> Pkg.add("QuantumOptics")
```

And will probably want to add Jupyter compatability by downloading IJulia:

```julia
julia> Pkg.add("IJulia")
```

after which you'll have the option to start a Julia kernel from Jupyter notebook.

### Updating IonSim

IonSim.jl is a work in progress. You can update your local code to the most recent version
with the following:

```julia
julia> using Pkg
julia> Pkg.update("IonSim")
```


[license-url]: https://github.com/HaeffnerLab/IonSim.jl/blob/master/LICENSE.md
[license-badge]: https://img.shields.io/badge/License-MIT-green.svg

[codecov-url]: https://codecov.io/gh/HaeffnerLab/IonSim.jl
[codecov-badge]: https://codecov.io/gh/HaeffnerLab/IonSim.jl/branch/master/graph/badge.svg

[travis-url]: https://travis-ci.org/HaeffnerLab/IonSim.jl
[travis-badge]: https://travis-ci.org/HaeffnerLab/IonSim.jl.svg?branch=master

[twitter-url]: https://twitter.com/Berkeley_ions
[twitter-badge]: https://img.shields.io/twitter/follow/Berkeley_ions.svg?style=social&label=@Berkeley_ions

[logo-url]: https://github.com/HaeffnerLab/IonSim.jl/blob/media/smallest_logo.png?raw=true
