using AJD
using Documenter

DocMeta.setdocmeta!(AJD, :DocTestSetup, :(using AJD); recursive = true)

makedocs(;
    modules = [AJD],
    authors = "Mustafa Erdenay Gürol erdenay16@gmail.com",
                "Johan Peter Fagrell peter.fagell@hotmail.se",
                "Jorinde Wellmer wellmer@campus.tu-berlin.de",
    sitename = "AJD.jl",
    format = Documenter.HTML(;
        canonical = "https://erdenay16.github.io/AJD.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/erdenay16/AJD.jl", devbranch = "main")
