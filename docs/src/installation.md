# Installation

If you haven't done so already, download [Julia](https://julialang.org/) (platform specific instructions can be found [here](https://julialang.org/downloads/)). Next, open the Julia app, which should launch a [REPL](https://docs.julialang.org/en/v1/stdlib/REPL/#The-Julia-REPL-1) session, and install IonSim using the following commands:

```julia
using Pkg
Pkg.add("IonSim")
```

The main way you'll want to interact with IonSim is inside of a [Jupyter notebook](https://jupyter.org/). This requires [IJulia.jl](https://github.com/JuliaLang/IJulia.jl):

```julia
Pkg.add("IJulia")
```

And finally, if you'd like to follow along with the example notebooks, you'll need at least [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl), an interface to Python's [Matplotlib](https://matplotlib.org/) library:

```julia
Pkg.add("PyPlot")
```

Updating to the latest version of IonSim is easy:

```julia
Pkg.update("IonSim")
```