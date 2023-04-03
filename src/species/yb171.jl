export Yb171

const properties_yb171 = SpeciesProperties(
    shortname="¹⁷¹Yb⁺",
    mass=28.8384644689030595108e-26,
    charge=1,
    nuclearspin=1//2,
    level_structure=OrderedDict(
        ls"6²S_1/2(F=0)"    => -9.48210909315e9,
        ls"6²S_1/2(F=1)"    => 3.16070303105e9,
        ls"5²D_3/2(F=0)"    => 6.88342460964640e14,
        ls"5²D_3/2(F=1)"    => 6.88350470564640e14,
        ls"5²D_5/2(F=2)"    => 7.29477895985202e14,
        ls"5²D_5/2(F=3)"    => 7.29474121985202e14,
        ls"6²P_1/2(F=0)"    => 8.112913749003559e14,
        ls"6²P_1/2(F=1)"    => 8.112934798003559e14,
        jk"²[3/2]_1/2(F=0)" => 1.008917341058788e15,
        jk"²[3/2]_1/2(F=1)" => 1.008917341058788e15,
        jk"²[3/2]_3/2(F=1)" => 8.621425511314839e14,
        jk"²[3/2]_3/2(F=2)" => 8.621425511314839e14,
        jk"²[5/2]_5/2(F=2)" => 9.70461163716380e14,
        jk"²[5/2]_5/2(F=3)" => 9.70461163716380e14,
        ls"6²F_7/2(F=3)"    => 6.42115934728750e14,
        ls"6²F_7/2(F=4)"    => 6.42119554728750e14,
    ),
    transitions=Dict(
        # Laser lines
        # 411nm
        (ls"6²S_1/2", ls"5²D_5/2")    => (multipole="E2", einsteinA=22),
        # 369nm
        (ls"6²S_1/2", ls"6²P_1/2")    => (multipole="E1", einsteinA=1.155e8),
        # 935nm
        (ls"5²D_3/2", jk"²[3/2]_1/2") => (multipole="E1", einsteinA=1.2e5),
        # 760nm
        (ls"6²F_7/2", jk"²[3/2]_3/2") => (multipole="E2", einsteinA=5e4),
        # 638nm
        # 976nm
        # 861nm
        # Decay lines
        # P1/2 -> D3/2
        (ls"6²P_1/2", ls"5²D_3/2")    => (multipole="E1", einsteinA=5.7e5),
        # [3/2]1/2 -> S1/2
        (ls"6²S_1/2", jk"²[3/2]_1/2") => (multipole="E1", einsteinA=8.05e7),
        # [3/2]3/2 -> S1/2
        (ls"6²P_1/2", jk"²[3/2]_3/2") => (multipole="E1", einsteinA=5.125e7),
    ),
    # Optional fields
    gfactors=Dict(
        ls"6²S_1/2(F=0)"    => 1.998,
        ls"6²S_1/2(F=1)"    => 1.998,
        ls"5²D_5/2(F=2)"    => 1.202,
        ls"5²D_5/2(F=3)"    => 1.202,
        ls"6²F_7/2(F=3)"    => 1.145,
        ls"6²F_7/2(F=4)"    => 1.145,
        jk"²[3/2]_1/2(F=0)" => 1.32,
        jk"²[3/2]_1/2(F=1)" => 1.32,
        jk"²[3/2]_3/2(F=1)" => 1.44,
        jk"²[3/2]_3/2(F=2)" => 1.44,
    ),
    nonlinear_zeeman=Dict(
        (ls"6²S_1/2(F=0)", 0) => B -> -155.305 * B^2,
        (ls"6²S_1/2(F=1)", 0) => B -> 155.305 * B^2
    )
)

# boilerplate code
IonInstance{:Yb171}(selected_sublevels=nothing, manualshift=Dict()) =
    IonInstance{:Yb171}(properties_yb171, selected_sublevels, manualshift)

Yb171 = IonInstance{:Yb171}
