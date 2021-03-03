push!(LOAD_PATH,"../src/")

using Documenter, SbmlToolkit

makedocs(
    sitename= "SbmlToolkit",    
    doctest = false,
    strict = false,
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md"
    ]
)
