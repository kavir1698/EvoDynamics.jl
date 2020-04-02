
# # Simplest model

# We can create and run simple Wright-Fisher simulations with EvoDynamics.jl. To that end, we define a single haploid species, in single region, with a single gene affecting a single phenotype. The values of parameters below are set arbitrarily.

using EvoDynamics

parameters = Dict(
  :L => (1),
  :P => (1),
  :R => (0.7),
  :C => nothing,
  :B => [[true]],
  :A =>  [[1.0]],
  :Q => [[1.0]],
  :Y => (0.5),
  :m => (1),
  :T => [[2.4]],
  :Ω => [[0.8]],
  :M => [(0.1, 0.0, 0.0)],
  :D => [(0.05, 0.0, 0.01)],
  :N => Dict(1 => (100)),
  :K => Dict(1 => [1000]),
  :migration_rates => [nothing],
  :E => (0.01),
  :generations => 10,
  :space => nothing
)

data, model = runmodel(parameters)
