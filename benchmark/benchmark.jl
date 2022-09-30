using EvoDynamics
using BenchmarkTools
using Agents

param_file = "benchmark/params.jl" 
param_file2 = "benchmark/params2.jl"

test_generations = [10, 20, 40, 80, 160, 320, 640, 1280]

results_time = Float64[]
results_memory = Float64[]
for gen in test_generations
  paramstring = readlines(param_file)
  paramstring[1] = "generations = $gen"
  open(param_file2, "w") do ff
    for ll in paramstring
      println(ff, ll)
    end
  end
  a = @benchmarkable run!(model, EvoDynamics.agent_step!, EvoDynamics.model_step!, $gen, adata=nothing, mdata=nothing, agents_first=false, showprogress=false) setup = (model = model_initiation(load_parameters(param_file2))) samples = 10 seconds = 120
  b = run(a)
  push!(results_time, minimum(b).time)
  push!(results_memory, minimum(b).memory)
end

#=
julia> print(results_time)
[6.891629e8, 5.5187392e9, 2.10339285e10, 2.56171767e10, 3.01109163e10, 2.81858753e10, 2.05297875e10, 3.37371201e10]
julia> print(results_memory)
[6.11760848e8, 4.15995912e9, 1.5822722048e10, 1.9975087205e10, 2.2654051072e10, 2.13194824e10, 1.5583138037e10, 2.4708061408e10]

using CairoMakie

fig = Figure(resolution = (800, 300))
ax = Axis(fig[1, 1], xlabel="Generations", ylabel="Time")
lines!(ax, test_generations, results_time)
ax2 = Axis(fig[2, 1], xlabel="Generations", ylabel="Memory")
lines!(ax2, test_generations, results_memory)
save("benchmark.pdf", fig)
=#