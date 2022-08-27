
# # Predator prey

# Here we create a predator prey model of two species, one haploid and one diploid.

using EvoDynamics

# Model parameters are in a .jl file as follows:

# ```julia
# ## 1. Functions

# function bn(agent::AbstractAgent, model::ABM)
#   return false
# end

# function optphens2(site::Tuple{Int,Int}, model::ABM)
#   opt1 = [0.8214794136627462, 0.49335401288317304, 0.3435556249656141, 0.25050033180075226]
#   opt2 = opt1 .+ 0.1
#   ss = EvoDynamics.coord2vertex(site, model)
#   if model.step[1] > 6
#     return [opt1[ss]]
#   else
#     return [opt2[ss]]
#   end
# end

# function optphens3(site::Tuple{Int,Int}, model::ABM)
#   phenotypic_values = [[0.7614758101208934 2.3911353663016586; 2.2091361414343313 1.7858661540288683; 0.7920974352892358 0.7630236717263839; 2.587205421882114 2.311211631439866],
#     [0.6437641305315445 2.097629262577973; 0.9954545836033715 1.1314248848669466; 2.469792530385348 1.490299526620522; 1.6158867433451882 0.2566056477862022]]
#   agent_site = (site[2] - 1) * size(model.space)[2] + (site[1])
#   if model.step[1] > 7
#     return phenotypic_values[1][agent_site, :]
#   else
#     return phenotypic_values[2][agent_site, :]
#   end
# end

# env_resources2(time::Int) = [200 220; 183 190]

# ## 2. Species parameters

# species1 = Dict(
#   :name => "a",
#   :number_of_genes => 7,
#   :number_of_phenotypes => 4,
#   :abiotic_phenotypes => [1],
#   :biotic_phenotypes => [2, 3],
#   :migration_phenotype => 4,
#   :migration_threshold => 3.5,
#   :vision_radius => 1,
#   :check_fraction => 0.5,
#   :ploidy => 2,
#   :epistasis_matrix => [1.0 0.43 -0.41 0.38 0.48 -0.43 -0.1 -0.08 -0.09 -0.5 0.41 0.44 -0.21 -0.12; 0.34 1.0 -0.19 -0.36 0.38 -0.28 0.24 -0.22 0.12 0.12 -0.12 -0.39 0.21 0.26; 0.05 0.27 1.0 0.04 0.01 -0.14 0.3 -0.28 0.43 -0.13 0.2 -0.02 0.25 -0.39; -0.12 0.33 -0.48 1.0 -0.4 -0.48 -0.22 -0.36 -0.24 -0.07 -0.12 -0.49 -0.37 0.27; 0.25 0.25 -0.14 0.49 1.0 0.28 -0.34 -0.49 0.45 -0.14 0.26 -0.13 -0.44 -0.17; -0.47 0.19 -0.24 0.41 0.08 1.0 0.11 0.03 0.15 0.49 0.04 0.41 -0.19 0.13; 0.37 0.09 -0.11 0.4 0.42 0.45 1.0 -0.01 -0.47 0.07 0.5 0.44 -0.18 -0.2; -0.32 0.15 0.4 -0.24 -0.21 0.5 0.22 1.0 -0.33 0.48 -0.49 0.07 0.5 -0.07; 0.02 -0.16 0.33 0.48 -0.42 0.39 0.2 -0.11 1.0 0.46 -0.06 0.22 -0.3 0.31; 0.41 -0.18 -0.16 -0.4 0.01 0.04 0.07 0.2 -0.37 1.0 -0.33 0.49 0.05 -0.42; 0.03 0.25 0.14 -0.36 0.28 -0.18 0.09 -0.2 0.46 -0.48 1.0 -0.21 0.41 -0.46; 0.0 0.44 -0.34 -0.42 0.37 -0.04 0.43 -0.25 0.21 0.19 0.29 1.0 -0.02 0.06; 0.46 -0.1 0.14 -0.22 -0.26 0.13 -0.5 -0.41 -0.31 0.0 -0.15 0.29 1.0 0.17; -0.22 0.21 0.46 -0.01 -0.35 -0.11 0.25 -0.03 0.18 -0.38 -0.4 -0.28 0.05 1.0],
#   :pleiotropy_matrix => Bool[1 1 0 0 0 0 0 1 0 1 0 1 0 1; 1 0 1 1 1 0 1 1 1 0 0 1 1 0; 1 0 1 0 1 0 0 0 0 0 1 0 1 0; 0 0 0 0 1 0 0 1 1 1 1 0 1 0],
#   :growth_rate => 0.8,
#   :expression_array => [0.28878032859775615, 0.4629421231828499, 0.26092147517051467, 0.952859489607121, 0.9638502824424, 0.05038142018016245, 0.05930756376654234, 0.033459292878885716, 0.32421526342800044, 0.9029235877297073, 0.7670060809312949, 0.12766808941531993, 0.8656895869985795, 0.342191940658253],
#   :selection_coefficient => 0.5,
#   :mutation_probabilities => [0.9, 0.0, 0.0],
#   :mutation_magnitudes => [0.05, 0.0, 0.01],
#   :N => [1000, 0, 0, 0],
#   :environmental_noise => 0.01,
#   :optimal_phenotypes => optphens2,
#   :bottleneck_function => bn,
#   :age => 5,
#   :recombination => 1,
#   :initial_energy => 0,
#   :reproduction_start_age => 1,
#   :reproduction_end_age => 5,
# )

# species2 = Dict(
#   :name => "b",
#   :number_of_genes => 8,
#   :number_of_phenotypes => 5,
#   :abiotic_phenotypes => [1, 2],
#   :biotic_phenotypes => [3, 4],
#   :migration_phenotype => 5,
#   :migration_threshold => 3.4,
#   :vision_radius => 1,
#   :check_fraction => 0.5,
#   :ploidy => 1,
#   :epistasis_matrix => [1.0 -0.01 -0.01 -0.34 -0.42 0.18 -0.09 0.13; 0.28 1.0 -0.21 -0.11 0.12 -0.46 0.21 -0.39; 0.0 0.28 1.0 0.1 0.09 0.3 0.1 0.17; 0.0 -0.42 0.05 1.0 -0.11 -0.07 -0.27 0.43; 0.31 0.47 -0.42 -0.34 1.0 -0.4 0.16 0.11; -0.19 0.29 -0.33 -0.49 0.17 1.0 -0.1 -0.28; -0.43 0.38 -0.06 0.39 0.21 -0.5 1.0 -0.08; 0.26 0.27 -0.44 -0.08 0.47 -0.27 -0.27 1.0],
#   :pleiotropy_matrix => Bool[1 0 0 0 1 1 0 0; 1 0 0 1 1 0 1 1; 0 1 1 1 1 1 0 1; 1 0 1 0 0 1 0 1; 1 1 0 1 1 1 1 1],
#   :growth_rate => 1.2,
#   :expression_array => [0.24923147816626035, 0.7155732641738595, 0.9655184311211502, 0.8638149724268844, 0.5075272565823061, 0.9189652626508431, 0.7897640036022151, 0.17091233765481717],
#   :selection_coefficient => 0.5,
#   :mutation_probabilities => [0.9, 0.0, 0.0],
#   :mutation_magnitudes => [0.05, 0.0, 0.01],
#   :N => [1000, 0, 0, 0],
#   :environmental_noise => 0.01,
#   :optimal_phenotypes => optphens3,
#   :bottleneck_function => bn,
#   :age => 3,
#   :recombination => 1,
#   :initial_energy => 0,
#   :reproduction_start_age => 1,
#   :reproduction_end_age => 3,
# )

# ## 3. Model parameters

# #NB this dict should be called model_parameters
# model_parameters = Dict(
#   :species => [species1, species2],
#   :generations => 14,
#   :space => (2, 2),
#   :metric => "chebyshev",
#   :periodic => false,
#   :resources => env_resources2,
#   :interactions => [0.1 0.0; 0.0 -1.0],
#   :food_sources => [1.0 0.5; 0.0 0.0],
#   :seed => nothing
# )
# ```

param_file = "../../examples/paramfile2.jl"  #hide
agentdata, modeldata, model = runmodel(param_file);

modeldata