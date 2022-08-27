
# # Simple Wright-Fisher

# We can create and run simple Wright-Fisher simulations with EvoDynamics.jl. To that end, we define a single haploid species, in an unstructured space, with two single genes affecting biotic and abiotic traits, respectively. 

using EvoDynamics

# A simple one-species model with no spatial structure. Model parameters are in a .jl file as follows:

# ```julia
# ## 1. Functions
# function bn(agent::AbstractAgent, model::ABM)
#   return false
# end

# function optphens(site::Tuple{Int,Int}, model::ABM)
#   return [1.5]
# end

# env_resources(time::Int) = [200]

# ## 2. Species parameters

# species1 = Dict(
#   :name => "a",
#   :number_of_genes => 2,
#   :number_of_phenotypes => 2,
#   :abiotic_phenotypes => [1],
#   :biotic_phenotypes => [2],
#   :migration_phenotype => 0,
#   :migration_threshold => 3.5,
#   :vision_radius => 0,
#   :check_fraction => 0,
#   :ploidy => 1,
#   :epistasis_matrix => [1.0 0.0; 0.0 1.0],
#   :pleiotropy_matrix => Bool[1 0; 0 1],
#   :growth_rate => 1.0,
#   :expression_array => [0.28, 0.46],
#   :selection_coefficient => 0.5,
#   :mutation_probabilities => [0.9, 0.0, 0.0],
#   :mutation_magnitudes => [0.05, 0.0, 0.0],
#   :N => [100],
#   :environmental_noise => 0.01,
#   :optimal_phenotypes => optphens,
#   :bottleneck_function => bn,
#   :age => 2,
#   :recombination => 0,
#   :initial_energy => 0
#   :reproduction_start_age => 1,
#   :reproduction_end_age => 1,
# )

# ## 3. Model parameters

# #NB this dict should be called model_parameters
# model_parameters = Dict(
#   :species => [species1],
#   :generations => 5,
#   :space => nothing,
#   :metric => "chebyshev",
#   :periodic => false,
#   :resources => env_resources,
#   :interactions => [-0.1],
#   :food_sources => [1.0],
#   :seed => nothing
# )
# ```

param_file = "../../examples/paramfile1.jl" #hide
agentdata, modeldata, model = runmodel(param_file);

modeldata