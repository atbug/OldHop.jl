using Documenter, Hop

makedocs()

deploydocs(
    repo = "github.com/mistguy/Hop.jl.git",
    deps = Deps.pip("mkdocs", "python-markdown-math"),
    julia = "1.0"
)
