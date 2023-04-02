export Ca40p

const properties_ca40 = 
    SpeciesProperties(;
        # shortname is a string identifier for the ion that gets used, e.g. during print 
        # statements. Any valid string will do
        shortname="⁴⁰Ca⁺",
        mass=6.635943757345042e-26,  # The mass of the ion in kg
        charge=1,                    # The net charge of the ion in units of fundamental charge
        nuclearspin=0,               # The total spin of the nucleus
        

        level_structure=OrderedDict(
            "S1/2" => (ls"4²S_1/2", 0),
            "D3/2" => (ls"3²D_3/2", 4.09335071228e14),
            "D5/2" => (ls"3²D_5/2", 4.1115503183857306e14),
            "P1/2" => (ls"4²P_1/2", 7.554e14),
            "P3/2" => (ls"4²P_3/2", 7.621e14),
        ),

        level_structure=OrderedDict(
            "S1/2" => (LSTerm(n=4, s=2, l=0, j=1//2), E=0),
            "D3/2" => (LSTerm(n=3, s=2, l=2, j=3//2), E=4.09335071228e14),
            "D5/2" => (LSTerm(n=3, s=2, l=2, j=5//2), E=4.1115503183857306e14),
            "P1/2" => (LSTerm(n=4, s=2, l=1, j=1//2), E=7.554e14),
            "P3/2" => (LSTerm(n=4, s=2, l=1, j=3//2), E=7.621e14),
        ),
        
        
        level_structure=OrderedDict(
            ls"4²S_1/2" => 0,
            ls"3²D_3/2" => 4.09335071228e14,
            ls"3²D_5/2" => 4.1115503183857306e14,
            ls"3²P_1/2" => 7.554e14,
            ls"3²P_3/2" => 7.621e14,
        ),

        full_transitions=Dict(
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