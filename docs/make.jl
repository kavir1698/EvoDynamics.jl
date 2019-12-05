using Pkg
Pkg.activate(@__DIR__)

using Documenter, Agents, Distributions, Random, StatsBase
using Literate

# %% Literate convertion
indir = joinpath(@__DIR__, "..", "examples")
outdir = joinpath(@__DIR__, "src")
for file in ("example1.jl", "example2.jl")
	Literate.markdown(joinpath(indir, file), outdir)
end

# %%
makedocs(
sitename= "EvoDynamics.jl",
authors = "Ali R. Vahdati and Carlos Melian.",
doctest = false,
format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    ),
pages = [
    "Introduction" => "index.md",
	"Tutorial" => "tutorial.md",
	"API" => "api.md",
	"Examples" => [
	  "Example 1" => "example1.md",
	  "Example 2" => "example2.md",
		],
    ],
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(repo = "github.com/kavir1698/EvoDynamics.git",
               target = "build")
end
