generations = 20
space = (1, 1)

## 1. Species parameters

species1 = Dict(
  :name => "a",
  :number_of_genes =>  2,
  :number_of_phenotypes => 2,
  :abiotic_phenotypes => [1],
  :biotic_phenotypes => [2],
  :migration_phenotype => 0,
  :migration_threshold => 3.5,
  :vision_radius => 0,
  :check_fraction => 0.0,
  :ploidy => 1,
  :epistasis_matrix => [1.0 0.0; 0.0 1.0],
  :pleiotropy_matrix => Bool[1 0; 0 1],
  :growth_rate => 1.0,
  :expression_array => [0.28, 0.46],
  :selection_coefficient => 0.1,
  :phenotype_contribution_to_fitness => nothing,
  :mutation_probabilities => [0.9, 0.0, 0.0],
  :mutation_magnitudes => [0.05, 0.0, 0.0],
  :N => [100],
  :environmental_noise => 0.01,
  :optimal_phenotypes => [fill([1.5 for p in 1:1], space...) for t in 0:generations],
  :bottlenecks => [fill(0.0, space...) for t in 0:generations],
  :age => 2,
  :recombination => 0,
  :initial_energy => 0,
  :reproduction_start_age => 1,
  :reproduction_end_age => 2,
  :abiotic_variance => 1.0,
  :biotic_variance => 1.0,
  :mating_scheme => 1
)

## 2. Model parameters

#NB this dict should be called model_parameters
model_parameters = Dict(
  :species => [species1],
  :generations => generations,
  :space => space,
  :metric => "chebyshev",
  :periodic => false ,
  :resources => [fill(200.0, 1,1) for i in 0:generations],
  :interactions => [-0.1],
  :food_sources => [1.0],
  :seed => nothing
)