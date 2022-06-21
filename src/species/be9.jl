using .PhysicalConstants: PhysicalConstant

export Be9

const properties_be9 = IonProperties(
    shortname = "⁹Be",
    mass = 1.496508080073e-26,
    charge = 1,
    nuclearspin = 3 // 2,
    full_level_structure = OrderedDict(
        "S1/2f=1" => (n = 2, l = 0, j = 1 // 2, f = 1, E = 0.78126104631e9),
        "S1/2f=2" => (n = 2, l = 0, j = 1 // 2, f = 2, E = -0.468756627786e9),
        "P1/2f=1" => (n = 2, l = 1, j = 1 // 2, f = 1, E = 957.4772214787497e12),
        "P1/2f=2" => (n = 2, l = 1, j = 1 // 2, f = 2, E = 957.4769842787498e12),
        "P3/2f=0" => (n = 2, l = 1, j = 3 // 2, f = 0, E = 957.6770501980199e12),
        "P3/2f=1" => (n = 2, l = 1, j = 3 // 2, f = 1, E = 957.6770491780198e12),
        "P3/2f=2" => (n = 2, l = 1, j = 3 // 2, f = 2, E = 957.6770471380199e12),
        "P3/2f=3" => (n = 2, l = 1, j = 3 // 2, f = 3, E = 957.6770440780199e12),
    ),
    full_transitions = Dict(
        ("S1/2f=1", "P1/2f=1") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=1", "P1/2f=2") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=2", "P1/2f=1") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=2", "P1/2f=2") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=1", "P3/2f=0") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=1", "P3/2f=1") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=1", "P3/2f=2") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=2", "P1/2f=1") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=2", "P1/2f=2") => (multipole = "E1", einsteinA = 19.4e6),
        ("S1/2f=2", "P1/2f=3") => (multipole = "E1", einsteinA = 19.4e6),
    ),

    # Optional fields
    default_sublevel_selection = [
        ("S1/2f=1", "all"),
        ("S1/2f=2", "all"),
        ("P1/2f=1", "all"),
        ("P1/2f=2", "all"),
        ("P3/2f=0", "all"),
        ("P3/2f=1", "all"),
        ("P3/2f=2", "all"),
        ("P3/2f=3", "all"),
    ],
)

# boilerplate code
IonInstance{:Be9}(
    selected_sublevels::Union{Vector{Tuple{String, T}}, String, Nothing} where {T} = nothing,
    starkshift = Dict()
) = IonInstance{:Be9}(properties_be9, selected_sublevels, starkshift)

Be9 = IonInstance{:Be9}