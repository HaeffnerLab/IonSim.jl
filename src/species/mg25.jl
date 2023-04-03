export Mg25

const properties_mg25 = SpeciesProperties(
    shortname="²⁵Mg⁺",
    mass=4.1489954812e-26,
    charge=1,
    nuclearspin=5//2,
    level_structure=OrderedDict(
        ls"3²S_1/2(F=2)" => 1.043445158e9,
        ls"3²S_1/2(F=3)" => -0.745317970e9,
        ls"3²P_1/2(F=2)" => 1069.3408690519877e12,
        ls"3²P_1/2(F=3)" => 1069.3405639519879e12,
        ls"3²P_3/2(F=1)" => 1072.0853411194748e12,
        ls"3²P_3/2(F=2)" => 1072.0853033394746e12,
        ls"3²P_3/2(F=3)" => 1072.0852466694748e12,
        ls"3²P_3/2(F=4)" => 1072.0851711094747e12,
    ),
    transitions=Dict(
        (ls"2²S_1/2", ls"2²P_1/2") => (multipole="E1", einsteinA=41.3e6),
        (ls"2²S_1/2", ls"2²P_3/2") => (multipole="E1", einsteinA=41.3e6),
    ),
    # Optional fields
)

# boilerplate code
IonInstance{:Mg25}(selected_sublevels=nothing, manualshift=Dict()) =
    IonInstance{:Mg25}(properties_mg25, selected_sublevels, manualshift)

Mg25 = IonInstance{:Mg25}
