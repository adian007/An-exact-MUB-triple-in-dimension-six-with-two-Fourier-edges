using Documenter, SLEEFPirates

makedocs(;
  modules = [SLEEFPirates],
  format = Documenter.HTML(),
  pages = ["Home" => "index.md"],
  sitename = "SLEEFPirates.jl",
  authors = "Chris Elrod",
)

deploydocs(; repo = "github.com/JuliaSIMD/SLEEFPirates.jl")
