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
##     5D₅/₂  →  5674.824 cm⁻¹
##     6P₁/₂  →  20261.563 cm⁻¹
##     6P₃/₂  →  21952.408 cm⁻¹
##   Converted to Hz via E(Hz) = E(cm⁻¹) × c × 100
##
## Einstein A coefficients:
##   E1 rates from Iskrenova-Tchoukova & Safronova, Phys. Rev. A 78, 012508 (2008)
##   E2 rates estimated from measured state lifetimes:
##     τ(5D₅/₂) ≈ 32 s  (Gudjons et al., Phys. Lett. A 220, 269 (1996))
##     τ(5D₃/₂) ≈ 80 s  (de Munshi et al., Phys. Rev. A 91, 040501 (2015))
##
## g-factors:
##   g(6S₁/₂) = 2.0024905  — theoretical (Ginges & Flambaum)
##   g(5D₅/₂) = 1.200371   — measured: Hoffman et al., Phys. Rev. A 88, 025401 (2013)
##   g(5D₃/₂) = 0.7993     — close to Landé value 4/5
##   g(6P₁/₂) = 0.6667     — Landé value 2/3
##   g(6P₃/₂) = 1.3333     — Landé value 4/3

const properties_ba138 = IonProperties(;
    shortname = "¹³⁸Ba",
    mass = 2.2899757e-25,       # 137.905247 u × 1.66053906660e-27 kg/u (AME2020)
    charge = 1,
    nuclearspin = 0,            # ¹³⁸Ba: I = 0, no hyperfine structure → f = j
    full_level_structure = OrderedDict(
        "S1/2" => (n=6, l=0, j=1 // 2, f=1 // 2, E=0),
        "D3/2" => (n=5, l=2, j=3 // 2, f=3 // 2, E=1.4611441e14),
        "D5/2" => (n=5, l=2, j=5 // 2, f=5 // 2, E=1.7012693e14),
        "P1/2" => (n=6, l=1, j=1 // 2, f=1 // 2, E=6.0742615e14),
        "P3/2" => (n=6, l=1, j=3 // 2, f=3 // 2, E=6.5857224e14),
    ),
    full_transitions = Dict(
        ("S1/2", "D5/2") => (multipole="E2", einsteinA=3.1e-2),   # 1762 nm; A ≈ 1/τ(D5/2)
        ("S1/2", "P1/2") => (multipole="E1", einsteinA=9.53e7),   # 493 nm  Doppler cooling
        ("D3/2", "P1/2") => (multipole="E1", einsteinA=3.10e7),   # 650 nm  repump
        ("S1/2", "D3/2") => (multipole="E2", einsteinA=1.2e-2),   # 2051 nm
        ("S1/2", "P3/2") => (multipole="E1", einsteinA=1.19e8),   # 455 nm
        ("D3/2", "P3/2") => (multipole="E1", einsteinA=6.0e6),
        ("D5/2", "P3/2") => (multipole="E1", einsteinA=4.12e7),   # 614 nm  shelving return
    ),
    default_sublevel_selection = [("S1/2", "all"), ("D5/2", "all")],
    gfactors = Dict(
        "S1/2" => 2.0024905,    # theoretical; close to free-electron g ≈ 2.00232
        "D5/2" => 1.200371,     # Hoffman et al., Phys. Rev. A 88, 025401 (2013)
        "D3/2" => 0.7993,       # Landé approximation (4/5 = 0.8)
        "P1/2" => 0.6667,       # Landé (2/3)
        "P3/2" => 1.3333,       # Landé (4/3)
    )
)

# boilerplate — identical pattern to ca40.jl, mg25.jl, yb171.jl
IonInstance{:Ba138}(selected_sublevels=nothing, manualshift=Dict()) =
    IonInstance{:Ba138}(properties_ba138, selected_sublevels, manualshift)
Ba138 = IonInstance{:Ba138}