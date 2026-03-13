using ThermoelasticAcoustic
using Documenter

DocMeta.setdocmeta!(
    ThermoelasticAcoustic, :DocTestSetup, :(using ThermoelasticAcoustic); recursive = true)

makedocs(;
    modules = [ThermoelasticAcoustic],
    authors = "Bruno Alves do Carmo <bruno.carmo@ppgi.ufrj.br>",
    sitename = "ThermoelasticAcoustic.jl",
    format = Documenter.HTML(;
        canonical = "https://bacarmo.github.io/ThermoelasticAcoustic.jl",
        edit_link = "main",
        assets = String[],
        size_threshold = 300 * 1024,  # error threshold: 300 KiB
        size_threshold_warn = 200 * 1024   # warning threshold: 200 KiB
    ),
    pages = [
        "Home" => "index.md",
        "Model" => "model.md",
        "Approximation Problem" => [
            "methods/scheme1.md",
            "methods/scheme2.md"
        ],
        "API" => "api.md"
    ]
)

deploydocs(;
    repo = "github.com/bacarmo/ThermoelasticAcoustic.jl",
    devbranch = "main"
)
