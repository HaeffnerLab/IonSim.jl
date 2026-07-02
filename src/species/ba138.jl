export Ba138

## ¹³⁸Ba⁺ — singly ionized barium-138
##
## ¹³⁸Ba is the most abundant stable isotope (71.7%) and has zero nuclear
## spin (bosonic), making it the standard choice for trapped-ion experiments.
## The 6S₁/₂ ↔ 5D₅/₂ E2 transition at 1762 nm serves as the optical qubit,
## analogous to the 4S₁/₂ ↔ 3D₅/₂ transition in ⁴⁰Ca⁺ at 729 nm.
##
## Level structure:
##   6S₁/₂  — ground state (qubit |0⟩)
##   5D₃/₂  — metastable, τ ≈ 80 s; repumped via 6P₁/₂ with 650 nm light
##   5D₅/₂  — metastable, τ ≈ 32 s (qubit |1⟩); addressed at 1762 nm
##   6P₁/₂  — dipole-allowed; Doppler cooling at 493 nm
##   6P₃/₂  — dipole-allowed; addressed at 455 nm
##
## Energy levels:
##   Wavenumbers from NIST Atomic Spectra Database (Ba II):
##     6S₁/₂  →  0        cm⁻¹
##     5D₃/₂  →  4873.852 cm⁻¹
##     5D₅/₂  →  5674.807 cm⁻¹
##     6P₁/₂  →  20261.561 cm⁻¹
##     6P₃/₂  →  21952.404 cm⁻¹
##   Converted to Hz via E(Hz) = E(cm⁻¹) × c × 100
##
## Einstein A coefficients:
##   E1 rates from Iskrenova-Tchoukova & Safronova, Phys. Rev. A 78, 012508 (2008)
##   E2 rates inferred from measured state lifetimes:
##     τ(5D₅/₂) ≈ 31.2(0.9)s  (Auchter et al., Phys. Rev. A 90, 060501(R) (2014))
##     τ(5D₃/₂) ≈ 79.8(4.6)s  (Yu, Nagourney, Dehmelt, Phys. Rev. Lett. 78, 4898 (1997))
##
## g-factors:
##   g(6S₁/₂) = 2.0024922(10)  — measured: Knöll et al., Phys. Rev. A 54, 1199 (1996)
##   g(5D₅/₂) = 1.200372(4)(7)   — measured: Hoffman et al., Phys. Rev. A 88, 025401 (2013)
##   g(5D₃/₂) = 0.7993278(3)     — measured: Knöll et al., Phys. Rev. A 54, 1199 (1996)
##   g(6P₁/₂) = 0.6667     — Landé value 2/3
##   g(6P₃/₂) = 1.3333     — Landé value 4/3

const properties_ba138 = IonProperties(;
    shortname = "¹³⁸Ba",
    mass = 2.28997051e-25,       # 137.905247 u × 1.66053906660e-27 kg/u (AME2020)
    charge = 1,
    nuclearspin = 0,            # ¹³⁸Ba: I = 0, no hyperfine structure → f = j
    full_level_structure = OrderedDict(
        "S1/2" => (n=6, l=0, j=1 // 2, f=1 // 2, E=0),
        "D3/2" => (n=5, l=2, j=3 // 2, f=3 // 2, E=1.46114407e14),
        "D5/2" => (n=5, l=2, j=5 // 2, f=5 // 2, E=1.70126434e14),
        "P1/2" => (n=6, l=1, j=1 // 2, f=1 // 2, E=6.07426318e14),
        "P3/2" => (n=6, l=1, j=3 // 2, f=3 // 2, E=6.58116515e14),
    ),
    full_transitions = Dict(
        ("S1/2", "D5/2") => (multipole="E2", einsteinA=3.21e-2),   # 1762 nm; A ≈ 1/τ(D5/2) 
        # τ = 31.2(9) s, Auchter et al., Phys. Rev. A 90, 060501(R) (2014)
        ("S1/2", "P1/2") => (multipole="E1", einsteinA=9.32e7),   # 493 nm  Doppler cooling
        ("D3/2", "P1/2") => (multipole="E1", einsteinA=3.41e7),   # 650 nm  repump
        ("S1/2", "D3/2") => (multipole="E2", einsteinA=1.25e-2),   # 2051 nm
        # τ = 79.8(4.6) s, Yu, Nagourney, Dehmelt, Phys. Rev. Lett. 78, 4898 (1997)
        ("S1/2", "P3/2") => (multipole="E1", einsteinA=1.18e8),   # 455 nm
        ("D3/2", "P3/2") => (multipole="E1", einsteinA=4.44e6),
        ("D5/2", "P3/2") => (multipole="E1", einsteinA=3.65e7),   # 614 nm  shelving return
    ),
    # P1/2, P3/2 branching fractions from:
    #   Arnold et al., "Measurements of the branching ratios for 6P1/2 decays in 138Ba+" (2019)
    #   Kurz et al., Phys. Rev. A 101, 062515 (2020) — P3/2 branching fractions
    default_sublevel_selection = [("S1/2", "all"), ("D5/2", "all")],
    gfactors = Dict(
        "S1/2" => 2.0024922,    # Knöll
        "D5/2" => 1.200372,     # Hoffman et al., Phys. Rev. A 88, 025401 (2013)
        "D3/2" => 0.7993278,    # Knöll et al., Phys. Rev. A 54, 1199 (1996)
        "P1/2" => 0.6667,       # Landé (2/3)
        "P3/2" => 1.3333,       # Landé (4/3)
    )
)

# boilerplate — identical pattern to ca40.jl, mg25.jl, yb171.jl
IonInstance{:Ba138}(selected_sublevels=nothing, manualshift=Dict()) =
    IonInstance{:Ba138}(properties_ba138, selected_sublevels, manualshift)
Ba138 = IonInstance{:Ba138}