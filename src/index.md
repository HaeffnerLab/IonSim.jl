# Introduction

```@raw html
<!-- Place this tag in your head or just before your close body tag. -->
<script async defer src="https://buttons.github.io/buttons.js"></script>
<!-- Place this tag where you want the button to render. -->
<a class="github-button" href="https://github.com/HaeffnerLab/IonSim.jl" data-color-scheme="no-preference: light_high_contrast; light: light_high_contrast; dark: dark_high_contrast;" data-icon="octicon-star" data-size="large" data-show-count="true" aria-label="Star HaeffnerLab/IonSim.jl on GitHub">Star</a>
```

IonSim.jl is a tool to simulate the dynamics of a collection of trapped ions interacting with laser light.

IonSim primarily performs two jobs:
1. Keep track of the physical parameters necessary for describing the system, with a structure and nomenclature designed to be intuitive for experimentalists.
2. Using these parameters, construct a function that quickly evaluates the system's Hamiltonian at a particular point in time.

This functional form of the Hamiltonian can then be used either as input to any of the solvers implemented in QuantumOptics.timeevolution, or in the native solver. The native solver is a thin wrapper around QuantumOptics functions that implement additional checks.

## Why use it?
IonSim is designed from an experimentalist's perspective. If you are looking to do efficient simulations of trapped ion systems, with high control of energy levels and vibrational modes, this is the package for you.

The code is MIT-licensed; find the code [on GitHub](https://github.com/HaeffnerLab/IonSim.jl).

something3
