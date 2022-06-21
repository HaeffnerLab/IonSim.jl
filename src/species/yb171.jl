using .PhysicalConstants: PhysicalConstant

export Yb171

const properties_yb171 = IonProperties(
    shortname = "¹⁷¹Yb",
    mass = 28.8384644689030595108e-26,
    charge = 1,
    nuclearspin = 1 // 2,
    full_level_structure = OrderedDict(
        # should this be negative, or should I increase the energy of every state by 9GHz?
        "S1/2f=0" => (n = 6, l = 0, j = 1 // 2, f = 0, E = -9.48210909315e9),
        "S1/2f=1" => (n = 6, l = 0, j = 1 // 2, f = 1, E = 3.16070303105e9),
        "D3/2f=0" => (n = 5, l = 2, j = 3 // 2, f = 1, E = 6.88342460964640e14),
        "D3/2f=1" => (n = 5, l = 2, j = 3 // 2, f = 2, E = 6.88350470564640e14),
        "D5/2f=2" => (n = 5, l = 2, j = 5 // 2, f = 2, E = 7.29477895985202e14),
        "D5/2f=3" => (n = 5, l = 2, j = 5 // 2, f = 3, E = 7.29474121985202e14),
        "P1/2f=0" => (n = 6, l = 1, j = 1 // 2, f = 0, E = 8.112913749003559e14),
        "P1/2f=1" => (n = 6, l = 1, j = 1 // 2, f = 1, E = 8.112934798003559e14),
        "[3/2]1/2f=0" =>
            (n = 6, l = nothing, j = 1 // 2, f = 0, E = 1.008917341058788e15),
        "[3/2]1/2f=1" =>
            (n = 6, l = nothing, j = 1 // 2, f = 1, E = 1.008917341058788e15),
        "[3/2]3/2f=1" =>
            (n = 6, l = nothing, j = 3 // 2, f = 1, E = 8.621425511314839e14),
        "[3/2]3/2f=2" =>
            (n = 6, l = nothing, j = 3 // 2, f = 2, E = 8.621425511314839e14),
        "[5/2]5/2f=2" =>
            (n = 6, l = nothing, j = 5 // 2, f = 2, E = 9.70461163716380e14),
        "[5/2]5/2f=3" =>
            (n = 6, l = nothing, j = 5 // 2, f = 3, E = 9.70461163716380e14),
        "F7/2f=3" => (n = 6, l = 3, j = 7 // 2, f = 3, E = 6.42115934728750e14),
        "F7/2f=4" => (n = 6, l = 3, j = 7 // 2, f = 4, E = 6.42119554728750e14),
    ),
    full_transitions = Dict(
        # Laser lines
        # 411nm
        ("S1/2f=0", "D5/2f=2") => (multipole = "E2", einsteinA = 22),
        ("S1/2f=1", "D5/2f=2") => (multipole = "E2", einsteinA = 22),
        ("S1/2f=1", "D5/2f=3") => (multipole = "E2", einsteinA = 22),
        # 369nm
        ("S1/2f=0", "P1/2f=1") => (multipole = "E1", einsteinA = 1.155e8),
        ("S1/2f=1", "P1/2f=0") => (multipole = "E1", einsteinA = 1.155e8),
        ("S1/2f=1", "P1/2f=1") => (multipole = "E1", einsteinA = 1.155e8),
        # 935nm
        ("D3/2f=1", "[3/2]1/2f=0") => (multipole = "E1", einsteinA = 1.2e5),
        ("D3/2f=1", "[3/2]1/2f=1") => (multipole = "E1", einsteinA = 1.2e5),
        ("D3/2f=2", "[3/2]1/2f=1") => (multipole = "E1", einsteinA = 1.2e5),
        # 760nm
        ("F7/2f=3", "[3/2]3/2f=1") => (multipole = "E2", einsteinA = 5e4),
        ("F7/2f=3", "[3/2]3/2f=2") => (multipole = "E2", einsteinA = 5e4),
        ("F7/2f=4", "[3/2]3/2f=2") => (multipole = "E2", einsteinA = 5e4),
        # 638nm
        # 976nm
        # 861nm
        # Decay lines
        # P1/2 -> D3/2
        ("P1/2f=0", "D3/2f=1") => (multipole = "E1", einsteinA = 5.7e5),
        ("P1/2f=1", "D3/2f=2") => (multipole = "E1", einsteinA = 5.7e5),
        ("P1/2f=1", "D3/2f=2") => (multipole = "E1", einsteinA = 5.7e5),
        # [3/2]1/2 -> S1/2
        ("[3/2]1/2f=0", "S1/2f=1") => (multipole = "E1", einsteinA = 8.05e7),
        ("[3/2]1/2f=1", "S1/2f=0") => (multipole = "E1", einsteinA = 8.05e7),
        ("[3/2]1/2f=1", "S1/2f=1") => (multipole = "E1", einsteinA = 8.05e7),
        # [3/2]3/2 -> S1/2
        ("[3/2]1/2f=0", "S1/2f=1") => (multipole = "E1", einsteinA = 5.125e7),
        ("[3/2]1/2f=1", "S1/2f=0") => (multipole = "E1", einsteinA = 5.125e7),
        ("[3/2]1/2f=1", "S1/2f=1") => (multipole = "E1", einsteinA = 5.125e7),
    ),
    # Optional fields
    default_sublevel_selection = [
        ("S1/2f=0", "all"),
        ("S1/2f=1", "all"),
        ("P1/2f=0", "all"),
        ("P1/2f=1", "all")
    ],
    gfactors = Dict(
        "S1/2f=0" => 1.998,
        "S1/2f=1" => 1.998,
        "D5/2f=2" => 1.202,
        "D5/2f=3" => 1.202,
        "F7/2f=3" => 1.145,
        "F7/2f=4" => 1.145,
        "[3/2]1/2f=0" => 1.32,
        "[3/2]1/2f=1" => 1.32,
        "[3/2]3/2f=1" => 1.44,
        "[3/2]3/2f=2" => 1.44,
    ),
    nonlinear_zeeman = Dict(
        ("S1/2f=0", 0) => B -> -155.305 * B^2,
        ("S1/2f=1", 0) => B -> 155.305 * B^2
    )
)

# boilerplate code
IonInstance{:Yb171}(
    selected_sublevels::Union{Vector{Tuple{String, T}}, String, Nothing} where {T} = nothing,
    starkshift = Dict()
) = IonInstance{:Yb171}(properties_yb171, selected_sublevels, starkshift)

Yb171 = IonInstance{:Yb171}