export Be9

const properties_be9 = SpeciesProperties(
    # scalar properties
    shortname="⁹Be⁺",
    mass=1.496508080073e-26,
    charge=1,
    nuclearspin=3 // 2,
    level_structure=OrderedDict(
        # S_1/2
        ls"2²S_1/2(F=1)" =>  0.78126104631e9,
        ls"2²S_1/2(F=2)" => -0.468756627786e9,
        # P_1/2
        ls"2²P_1/2(F=1)" => 957.2010729076436e12,
        ls"2²P_1/2(F=2)" => 957.2008357076436e12,
        # P_3/2
        ls"2²P_3/2(F=0)" => 957.3965669467407e12,
        ls"2²P_3/2(F=1)" => 957.3965659267407e12,
        ls"2²P_3/2(F=2)" => 957.3965638867406e12,
        ls"2²P_3/2(F=3)" => 957.3965608267406e12,
    ),

    transitions=Dict(
        (ls"2²S_1/2", ls"2²P_1/2") => (multipole="E1", einsteinA=19.4e6),
        (ls"2²S_1/2", ls"2²P_3/2") => (multipole="E1", einsteinA=19.4e6),
    ),
    # # Optional fields
    
    # gfactors=Dict("S1/2" => 2.00225664, "D5/2" => 1.2003340)

    # quadratic_zeeman = Dict(("S1/2", -1//2) => B->1.3e-4*B^2,
    #                        ("D5/2", -5//2) => B->4.5e-4*B^2)  # Syntax example, not numerically accurate
)

# boilerplate code
IonInstance{:Be9}(selected_sublevels=nothing, manualshift=Dict()) =
    IonInstance{:Be9}(properties_be9, selected_sublevels, manualshift)

Be9 = IonInstance{:Be9}
