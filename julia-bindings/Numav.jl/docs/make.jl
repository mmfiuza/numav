# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

using Documenter, Numav, Markdown

function parse_md_file(filename)
    return Markdown.parse(read(joinpath(@__DIR__, "src", filename), String))
end

makedocs(
    sitename="Numav.jl",
    # modules = [Numav], # error if some docstring is forgotten to be showed
    meta = Dict(:CollapsedDocStrings => true),
    format = Documenter.HTML(
        assets = ["assets/custom.css"],
        prettyurls = false,
    ),
    pages = [
        "Home" => "index.md",
        "Modules" => [
            "FEM" => [
                "Helmholtz" => "fem-helmholtz.md",
            ],
        ],
        "Frequency-dependent physical quantities" => "fdpq.md",
        "Credits" => "credits.md",
        "License" => "license.md",
    ],
)
