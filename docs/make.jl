using Documenter, IonSim
push!(LOAD_PATH, "/Users/josephbroz/Desktop/IonSim/IonSim.jl/src") 

builddir = "build"

pages = [
    "index.md",
    "installation.md",
    "tutorial.md",
    "Constructing a Trapped Ion Experiment" => [
        "objects/ions.md",
        "IonTraps" => [
            "objects/iontraps.md",
            "objects/linearchains.md",
            "objects/vibrationalmodes.md"
        ],
        "objects/chambers.md",
        "objects/hamiltonian.md",
    ],
    "timeevolution.md",
    "examples.md",
    "api.md"
]

makedocs(
    modules = [IonSim],
    checkdocs = :exports,
    format=Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = nothing,
        sidebar_sitename = false,
        canonical = "https://docs.ionsim.org/",
        assets = [
            asset("assets/js/link-back.js", class=:js, islocal=true),
            asset("assets/small-logo.png", class=:ico, islocal=true)
        ],
    ),
    build = builddir,
    sitename = "IonSim.jl",
    pages = pages
)

deploydocs(
    repo = "github.com/HaeffnerLab/IonSim.jl-Documentation.git",
    devurl = "",
    versions = nothing
)
