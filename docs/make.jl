#
# To create the html documentation locally, use
# the following command:
# 
#    julia --color=yes --project make.jl --local-html
#
# For the Latex pdf documentation, use instead:
#
#    julia --color=yes --project make.jl --latex-pdf
#

push!(LOAD_PATH,"../src/")

using Documenter
using ConleyDynamics
using DocumenterCitations
using Bibliography
using Pkg

cdv = Pkg.TOML.parsefile(joinpath(pkgdir(ConleyDynamics), "Project.toml"))["version"]

#bib = CitationBibliography(
#    joinpath(@__DIR__, "src", "refs.bib");
#    style=:alpha) # The default is     style=:numeric
#                  # One can also use   style=:authoryear
#
# sort_bibliography!(bib.entries, :nyt)  # name-year-title

pageslist = ["Overview" => "index.md",
             "Manual" => Any[
                 "man/tutorial.md",
                 "man/lefschetz.md",
                 "man/homology.md",
                 "man/conley.md",
                 "man/examples.md",
                 "man/sparse.md",
                 "man/references.md"
                  ],
             "Core API" => Any[
                 "apicore/datastruct.md",
                 "apicore/lefschetz.md",
                 "apicore/homology.md",
                 "apicore/conley.md",
                 "apicore/examples.md",
                 "apicore/plots.md",
                 "apicore/sparse.md",
                 "apicore/apiindex.md"
                  ],
            ]

DocMeta.setdocmeta!(ConleyDynamics, :DocTestSetup,
                    :(using ConleyDynamics); recursive=true)

if "--local-html" in ARGS
     makedocs(sitename="ConleyDynamics.jl",
        modules=[ConleyDynamics],
        format = Documenter.HTML(prettyurls = false,
                                 assets=String["assets/citations.css","assets/favicon.ico"]),
        pages = pageslist,
        authors = "Thomas Wanner",
        plugins = [CitationBibliography(
                     joinpath(@__DIR__, "src", "refs.bib");
                     style=:alpha) # The default is     style=:numeric
                                   # One can also use   style=:authoryear
                  ]
        )
elseif "--latex-pdf" in ARGS
    makedocs(sitename="ConleyDynamics.jl",
        modules=[ConleyDynamics],
        format = Documenter.LaTeX(platform = "docker", version=cdv),
        build = "latex_build",
        pages = pageslist,
        authors = "Thomas Wanner",
        plugins = [CitationBibliography(
                     joinpath(@__DIR__, "src", "refs.bib");
                     style=:alpha) # The default is     style=:numeric
                                   # One can also use   style=:authoryear
                  ]
        )
    cp(joinpath(@__DIR__, "latex_build", "ConleyDynamics.jl-$(cdv).pdf"),
        joinpath(@__DIR__, "build", "ConleyDynamics.pdf");
        force = true,
        )
else
    makedocs(sitename="ConleyDynamics.jl",
        modules=[ConleyDynamics],
        format = Documenter.HTML(assets=String["assets/citations.css","assets/favicon.ico"]),
        build = "build",
        pages = pageslist,
        authors = "Thomas Wanner",
        plugins = [CitationBibliography(
                     joinpath(@__DIR__, "src", "refs.bib");
                     style=:alpha) # The default is     style=:numeric
                                   # One can also use   style=:authoryear
                  ]
        )
    makedocs(sitename="ConleyDynamics.jl",
        modules=[ConleyDynamics],
        format = Documenter.LaTeX(platform = "docker", version=cdv),
        build = "latex_build",
        pages = pageslist,
        authors = "Thomas Wanner",
        plugins = [CitationBibliography(
                     joinpath(@__DIR__, "src", "refs.bib");
                     style=:alpha) # The default is     style=:numeric
                                   # One can also use   style=:authoryear
                  ]
        )
    cp(joinpath(@__DIR__, "latex_build", "ConleyDynamics.jl-$(cdv).pdf"),
        joinpath(@__DIR__, "build", "ConleyDynamics.pdf");
        force = true,
        )
    deploydocs(
        repo = "github.com/almost6heads/ConleyDynamics.jl.git",
        )
end

