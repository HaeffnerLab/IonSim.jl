export Mg25

const properties_mg25 = IonProperties(
    shortname = "²⁵Mg",
    mass = 4.1489954812e-26,
    charge = 1,
    nuclearspin = 5 // 2,
    full_level_structure = OrderedDict(
        "S1/2f=2" => (n = 3, l = 0, j = 1 // 2, f = 2, E = 1.043445158e9),
        "S1/2f=3" => (n = 3, l = 0, j = 1 // 2, f = 3, E = -0.745317970e9),
        "P1/2f=2" => (n = 3, l = 1, j = 1 // 2, f = 2, E = 1069.3408690519877e12),
        "P1/2f=3" => (n = 3, l = 1, j = 1 // 2, f = 3, E = 1069.3405639519879e12),
        "P3/2f=1" => (n = 3, l = 1, j = 3 // 2, f = 1, E = 1072.0853411194748e12),
        "P3/2f=2" => (n = 3, l = 1, j = 3 // 2, f = 2, E = 1072.0853033394746e12),
        "P3/2f=3" => (n = 3, l = 1, j = 3 // 2, f = 3, E = 1072.0852466694748e12),
        "P3/2f=4" => (n = 3, l = 1, j = 3 // 2, f = 4, E = 1072.0851711094747e12),
    ),
    full_transitions = Dict(
        ("S1/2f=2", "P1/2f=2") => (multipole = "E1", einsteinA = 41.3e6),
        ("S1/2f=2", "P1/2f=3") => (multipole = "E1", einsteinA = 41.3e6),
        ("S1/2f=3", "P1/2f=2") => (multipole = "E1", einsteinA = 41.3e6),
        ("S1/2f=3", "P1/2f=3") => (multipole = "E1", einsteinA = 41.3e6),
        ("S1/2f=2", "P3/2f=1") => (multipole = "E1", einsteinA = 41.3e6),
        ("S1/2f=2", "P3/2f=2") => (multipole = "E1", einsteinA = 41.3e6),
        ("S1/2f=2", "P3/2f=3") => (multipole = "E1", einsteinA = 41.3e6),
        ("S1/2f=3", "P1/2f=2") => (multipole = "E1", einsteinA = 41.3e6),
        ("S1/2f=3", "P1/2f=3") => (multipole = "E1", einsteinA = 41.3e6),
        ("S1/2f=3", "P1/2f=4") => (multipole = "E1", einsteinA = 41.3e6),
    ),
    # Optional fields
    default_sublevel_selection = [
        ("S1/2f=2", "all"),
        ("S1/2f=3", "all"),
        ("P1/2f=2", "all"),
        ("P1/2f=3", "all"),
        ("P3/2f=1", "all"),
        ("P3/2f=2", "all"),
        ("P3/2f=3", "all"),
        ("P3/2f=4", "all"),
    ],
)

# boilerplate code
IonInstance{:Mg25}(
    selected_sublevels::Union{Vector{Tuple{String, T}}, String, Nothing} where {T} = nothing,
    manualshift = Dict()
) = IonInstance{:Mg25}(properties_mg25, selected_sublevels, manualshift)

Mg25 = IonInstance{:Mg25}
