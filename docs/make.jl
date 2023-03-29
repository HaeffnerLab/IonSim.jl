using Documenter, IonSim

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
            "objects/vibrationalmodes.md",
        ],
        "objects/emfields.md",
        "objects/waveforms.md",
        "objects/chambers.md",
    ],
    "hamiltonian.md",
    "Simulation" => [
            "timeevolution/solve.md",
            "timeevolution/rotatingframes.md",
            "timeevolution/technicalnoise.md",   
            "timeevolution/quantumnoise.md",
    ],
    "examples.md",
    "api.md",
]

makedocs(
    sitename = "IonSim.jl",
    authors = "Joseph Broz",
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
    pages = pages
)

deploydocs(
    repo = "github.com/HaeffnerLab/IonSim.jl.git",
    push_preview = true
)
