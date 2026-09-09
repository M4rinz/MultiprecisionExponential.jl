using MultiprecisionExponential
using Documenter

DocMeta.setdocmeta!(MultiprecisionExponential, :DocTestSetup, :(using MultiprecisionExponential); recursive=true)

makedocs(;
    modules=[MultiprecisionExponential],
    authors="Andrea Marino <134485532+M4rinz@users.noreply.github.com> and contributors",
    sitename="MultiprecisionExponential.jl",
    format=Documenter.HTML(;
        canonical="https://M4rinz.github.io/MultiprecisionExponential.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/M4rinz/MultiprecisionExponential.jl",
    devbranch="main",
)
