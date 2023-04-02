export Ca40

const properties_ca40 = SpeciesProperties(;
    shortname="⁴⁰Ca⁺",
    mass=6.635943757345042e-26,
    charge=1,
    nuclearspin=0,
    level_structure=OrderedDict(
        "S1/2" => (n=4, l=0, j=1 // 2, f=1 // 2, E=0),
        "D3/2" => (n=3, l=2, j=3 // 2, f=3 // 2, E=4.09335071228e14),
        "D5/2" => (n=3, l=2, j=5 // 2, f=5 // 2, E=4.1115503183857306e14),
        "P1/2" => (n=4, l=1, j=1 // 2, f=1 // 2, E=7.554e14),
        "P3/2" => (n=4, l=1, j=3 // 2, f=3 // 2, E=7.621e14),
    ),
    transitions=Dict(
        ("S1/2", "D5/2") => (multipole="E2", einsteinA=8.562e-1),
        ("S1/2", "P1/2") => (multipole="E1", einsteinA=1.299e8),
        ("D3/2", "P1/2") => (multipole="E1", einsteinA=1.060e7),
        ("S1/2", "D3/2") => (multipole="E2", einsteinA=9.259e-1),
        ("S1/2", "P3/2") => (multipole="E1", einsteinA=1.351e8),
        ("D3/2", "P3/2") => (multipole="E1", einsteinA=1.110e6),
        ("D5/2", "P3/2") => (multipole="E1", einsteinA=9.901e6),
    ),

    # Optional fields
    gfactors=Dict("S1/2" => 2.00225664, "D5/2" => 1.2003340)

    #nonlinear_zeeman = Dict(("S1/2", -1//2) => B->1.3e-4*B^2,
    #                        ("D5/2", -5//2) => B->4.5e-4*B^2)  # Syntax example, not numerically accurate
)

# boilerplate code
IonInstance{:Ca40}(selected_sublevels=nothing, manualshift=Dict()) =
    IonInstance{:Ca40}(properties_ca40, selected_sublevels, manualshift)

Ca40 = IonInstance{:Ca40}