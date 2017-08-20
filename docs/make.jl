using Documenter, Hop

makedocs(
   format = :html,
   sitename = "Hop.jl",
   pages = [
        "index.md",
    ]
)

deploydocs(
    repo = "github.com/mistguy/Hop.jl.git"
)
