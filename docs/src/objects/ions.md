# Ion

The Ion is the simplest `struct` in IonSim.

We have pre-loaded parameters for several ions, including Be9, Ca40, Mg25, and Yb171 ([full list here](https://github.com/HaeffnerLab/IonSim.jl/tree/master/src/species)).

You can also configure the parameters for any ion you want.

We support
* E1 (electric dipole) transitions
* E2 (electric quadrapole) transitions
* (coming soon) M1 (magnetic dipole) transitions
* Zeeman shifts
* AC (dynamic) Stark shifts

We do not support DC Stark shifts.