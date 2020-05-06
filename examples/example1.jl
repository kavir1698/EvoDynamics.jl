
# # Simplest model

# We can create and run simple Wright-Fisher simulations with EvoDynamics.jl. To that end, we define a single haploid species, in single region, with a single gene affecting a single phenotype. The values of parameters below are set arbitrarily.

using EvoDynamics

parameters = Dict(
  :ngenes => (1),
  :nphenotypes => (1),
  :growthrates => (0.7),
  :competitionCoeffs => nothing,
  :pleiotropyMat => [[true]],
  :epistasisMat =>  [[1.0]],
  :expressionArrays => [[1.0]],
  :selectionCoeffs => (0.5),
  :ploidy => (1),
  :optPhenotypes => [[2.4]],
  :covMat => [[0.8]],
  :mutProbs => [(0.1, 0.0, 0.0)],
  :mutMagnitudes => [(0.05, 0.0, 0.01)],
  :N => Dict(1 => (100)),
  :K => Dict(1 => [1000]),
  :migration_rates => [nothing],
  :E => (0.01),
  :generations => 10,
  :space => nothing
)

agentdata, modeldata, model = runmodel(parameters)
