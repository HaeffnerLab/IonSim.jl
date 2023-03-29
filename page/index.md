<!-- =============================
     ABOUT
    ============================== -->

\begin{section}{title="", name="About"}


\lead{IonSim is built on top of [QuantumOptics.jl](https://qojulia.org/) and makes it easy to simulate the dynamics of trapped ion systems like ths one:

\figure{path="assets/trapped_ion_system.png"}

Described by the Hamiltonian:}

$$
$$

IonSim primarily performs two jobs:
1. Keeps track of the physical parameters necessary for describing the system, with a structure and nomenclature designed to be intuitive for experimentalists.
2. Using these parameters, constructs a function that quickly evaluates the system's Hamiltonian at a particular point in time.

This Hamiltonian can then be used either as input to any of the solvers implemented in QuantumOptics, or in the native solver, which supports a straightforward treatment of technical noise.
\end{section}

<!-- ==============================
     GETTING STARTED
     ============================== -->
\begin{section}{title="Getting started"}

Install IonSim with:

```julia-repl
julia> using Pkg
julia> Pkg.add("IonSim")
```

You'll also likely want to add:

```julia-repl
julia> Pkg.add("IJulia")
julia> Pkg.add("PyPlot")
```
The first line adds Jupyter compatibility, so you'll have the option to start a Julia kernel from Jupyter notebook. This is the suggested way to interact with IonSim.
The second line adds an interface to Python's [Matplotlib](https://matplotlib.org/) library, which will make it easy to follow along with the examples.

IonSim is a work in progress, but you can keep current with the latest updates easily with:

```julia-repl
julia> Pkg.update("IonSim")
```

Once you've installed IonSim you can take a look at the tutorial or the examples.

<!-- ==============================
     Examples
     ============================== -->
\begin{section}{title="Examples"}

<!-- \lead{
} -->

\end{section}


<!-- =============================
     Collaborate
    ============================== -->

\begin{section}{title="Collaborate"}

\lead{
    Franklin can run your Julia code on the fly and show the output.
}


\end{section}


<!-- =============================
     Cite
    ============================== -->

\begin{section}{title="Cite"}

\lead{Make your page available online easily by leveraging GitHub Actions and GitHub Pages.}

By following these instructions, the content of the rendered website will be copied to a `gh-pages` branch where it will be deployed by GitHub.
If you would like to deploy the page with your own URL or using something else than GitHub, have a look at the specific instructions further on.

**Adjust DeployPage**: start by checking the `.github/workflows/DeployPage.yml` in particular:
* if you want to use Python or matplotlib, uncomment the relevant lines
* in the `run` block ensure that
    * `NodeJS` and `PkgPage` are added,
    * any packages that your page might rely on are added,
    * the `optimize` call has the appropriate `input` and `output` path (if you're in the default setting, leave as is).

**GitIgnore**: it's important you specify that `page/__site` should be ignored by git and not pushed to your repository otherwise the build process might not work properly. To do so create a file `.gitignore` containing the line

```
page/__site
```

as shown [here](https://github.com/tlienart/PkgPage.jl/blob/cce098535eb95c2c3ba919d605792abfee57710c/.gitignore#L3).

**GitAttributes**: in order for GitHub to ignore `page` folder it the language statistics for your repository, make sure to add a file `.gitattributes` with content

```
page/* linguist-vendored
```

like [this](https://github.com/tlienart/PkgPage.jl/blob/master/.gitattributes).

Now whenever you push changes to the `master` branch of your package, the  build process will be triggered and your page updated and deployed.
**That's it**.

**Avoiding clashes with Documenter.jl**: if you already use [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) you might want your page to be deployed in a specific folder of `gh-pages` as Documenter also generates files in `gh-pages`.

\alert{This will typically not be necessary as the names created by PkgPage and Documenter don't clash, but you might still prefer to avoid mixing the two (in which case, read on).}

you can do so in two steps:

1. change the `run` part of `DeployPage.yml` by specifying the `output` keyword argument  in `PkgPage.optimize` for instance: `PkgPage.optimize(input="page", output="page")`,
1. change the `prepath` in `config.md` to reflect that the base URL will contain that additional folder, for instance `@def prepath = "YourPackage.jl/page"`.

**Use your own URL**: you can usually get host services like Netlify to deploy a specific branch of a GitHub repo, do make sure to set `@def prepath = ""` in your `config.md` though.

If you want to do the deployment without GitHub actions then you will need to:

* ensure you have `purgecss` and `highlights` installed and available to `NodeJS`, the simplest way to do this is to install them via `NodeJS` with

```
using NodeJS;
run(`$(npm_cmd()) install highlight.js`);
run(`$(npm_cmd()) install purgecss`);
```
\\
* run `PkgPage.optimize(input="page", output="")` (adapting `input` as required)
* place the content of `page/__site` wherever your server requires it.

\end{section}
