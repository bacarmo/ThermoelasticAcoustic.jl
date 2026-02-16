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
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Model" => "model.md",
        "Approximation Problem" => "methods/crank_nicolson_galerkin.md",
        "API" => "api.md"
    ]
)

deploydocs(;
    repo = "github.com/bacarmo/ThermoelasticAcoustic.jl",
    devbranch = "main"
)
