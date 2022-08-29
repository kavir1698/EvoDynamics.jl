generations = 14
space = (2, 2)

function env_resources(time::Int)
  if time < 7
    return [200 158; 183 190]
  else
    return [200 158; 183 190] .- 50
  end
end

## 1. Species properties
##----------------------------------------------------------------

species1 = Dict(
  :name => "a",
  :number_of_genes => 7,
  :number_of_phenotypes => 4,
  :abiotic_phenotypes => [1],
  :biotic_phenotypes => [2, 3],
  :migration_phenotype => 4,  # can be 0 for no migration
  :migration_threshold => 3.2,  # phenotypic threshold after which migration is possible
  :vision_radius => 1,  # the radius of neighboring sites that the agent can see
  :check_fraction => 0.5, # the fraction of the observable sites that the agent checks
  :ploidy => 2, # A random epistasis matrix where the diagonal is 1.0, meaning that each locus affects itself 100%.
  # (ngenes x ploidy)^2
  :epistasis_matrix => [1.0 0.43 -0.41 0.38 0.48 -0.43 -0.1 -0.08 -0.09 -0.5 0.41 0.44 -0.21 -0.12; 0.34 1.0 -0.19 -0.36 0.38 -0.28 0.24 -0.22 0.12 0.12 -0.12 -0.39 0.21 0.26; 0.05 0.27 1.0 0.04 0.01 -0.14 0.3 -0.28 0.43 -0.13 0.2 -0.02 0.25 -0.39; -0.12 0.33 -0.48 1.0 -0.4 -0.48 -0.22 -0.36 -0.24 -0.07 -0.12 -0.49 -0.37 0.27; 0.25 0.25 -0.14 0.49 1.0 0.28 -0.34 -0.49 0.45 -0.14 0.26 -0.13 -0.44 -0.17; -0.47 0.19 -0.24 0.41 0.08 1.0 0.11 0.03 0.15 0.49 0.04 0.41 -0.19 0.13; 0.37 0.09 -0.11 0.4 0.42 0.45 1.0 -0.01 -0.47 0.07 0.5 0.44 -0.18 -0.2; -0.32 0.15 0.4 -0.24 -0.21 0.5 0.22 1.0 -0.33 0.48 -0.49 0.07 0.5 -0.07; 0.02 -0.16 0.33 0.48 -0.42 0.39 0.2 -0.11 1.0 0.46 -0.06 0.22 -0.3 0.31; 0.41 -0.18 -0.16 -0.4 0.01 0.04 0.07 0.2 -0.37 1.0 -0.33 0.49 0.05 -0.42; 0.03 0.25 0.14 -0.36 0.28 -0.18 0.09 -0.2 0.46 -0.48 1.0 -0.21 0.41 -0.46; 0.0 0.44 -0.34 -0.42 0.37 -0.04 0.43 -0.25 0.21 0.19 0.29 1.0 -0.02 0.06; 0.46 -0.1 0.14 -0.22 -0.26 0.13 -0.5 -0.41 -0.31 0.0 -0.15 0.29 1.0 0.17; -0.22 0.21 0.46 -0.01 -0.35 -0.11 0.25 -0.03 0.18 -0.38 -0.4 -0.28 0.05 1.0],
  # nphenotype x (ngenes x ploidy)
  :pleiotropy_matrix => Bool[1 1 0 0 0 0 0 1 0 1 0 1 0 1; 1 0 1 1 1 0 1 1 1 0 0 1 1 0; 1 0 1 0 1 0 0 0 0 0 1 0 1 0; 0 0 0 0 1 0 0 1 1 1 1 0 1 0],
  :growth_rate => 0.8 , # max number of offsprings per mating mean of a Poisson
  :expression_array => [0.28878032859775615, 0.4629421231828499, 0.26092147517051467, 0.952859489607121, 0.9638502824424, 0.05038142018016245, 0.05930756376654234, 0.033459292878885716, 0.32421526342800044, 0.9029235877297073, 0.7670060809312949, 0.12766808941531993, 0.8656895869985795, 0.342191940658253],  # ngenes x ploidy long
  :selection_coefficient => 0.2,
  :mutation_probabilities => [0.9, 0.0, 0.0],  # for gene expression array, pleiotropy matrix and epistasis matrix, respectively
  :mutation_magnitudes => [0.05, 0.0, 0.01], # same as above
  :N => [1000, 0, 0, 0],  # number of individuals per site at time 0
  :environmental_noise => 0.0,  # variance of a normal distribution with mean 0
  # each row is the optimal phenotypes for each site for all abiotic traits. There are as many element as number of sites times number of abiotic traits. The first N elements are for the first abiotic trait, where N is the number of sites, and so on.
  :optimal_phenotypes => [fill([1.5 for p in 1:1], space...) for t in 0:generations], # optimal phenotypes per site and time for each phenotype
  :age => 5,  # max age
  :reproduction_start_age => 1,
  :reproduction_end_age => 5,
  :recombination => 1, # Mean of a Poisson distributions for number of crossing overs
  :initial_energy => 0, # A parameter for parental care of infants. Values more than 0 indicate that newly born individuals can survive for a number of times without requiring food from the environment/other species. The consumption rate (i.e. how many generations this initial energy suffices) is determined by the sum of the corresponding rows in "food sources"
  :bottleneck_function => [fill(0.0, space...) for t in 0:generations]  # an array of arrays with probablity of external death at each site and generation.
)

species2 = Dict(
  :name => "b",
  :number_of_genes => 8,
  :number_of_phenotypes => 5,
  :abiotic_phenotypes => [1, 2],
  :biotic_phenotypes => [3, 4],
  :migration_phenotype => 5,
  :migration_threshold => 3.4,
  :vision_radius => 1,
  :check_fraction => 0.5,
  :ploidy => 1,
  :epistasis_matrix => [1.0 -0.01 -0.01 -0.34 -0.42 0.18 -0.09 0.13; 0.28 1.0 -0.21 -0.11 0.12 -0.46 0.21 -0.39; 0.0 0.28 1.0 0.1 0.09 0.3 0.1 0.17; 0.0 -0.42 0.05 1.0 -0.11 -0.07 -0.27 0.43; 0.31 0.47 -0.42 -0.34 1.0 -0.4 0.16 0.11; -0.19 0.29 -0.33 -0.49 0.17 1.0 -0.1 -0.28; -0.43 0.38 -0.06 0.39 0.21 -0.5 1.0 -0.08; 0.26 0.27 -0.44 -0.08 0.47 -0.27 -0.27 1.0],
  :pleiotropy_matrix => Bool[1 0 0 0 1 1 0 0; 1 0 0 1 1 0 1 1; 0 1 1 1 1 1 0 1; 1 0 1 0 0 1 0 1; 1 1 0 1 1 1 1 1],
  :growth_rate => 1.2,
  :expression_array => [0.24923147816626035, 0.7155732641738595, 0.9655184311211502, 0.8638149724268844, 0.5075272565823061, 0.9189652626508431, 0.7897640036022151, 0.17091233765481717],
  :selection_coefficient => 0.5,
  :mutation_probabilities => [0.9, 0.0, 0.0],
  :mutation_magnitudes => [0.05, 0.0, 0.01],
  :N => [1000, 0, 0, 0],
  :environmental_noise => 0.01,
  :optimal_phenotypes => [fill([1.5 for i in 1:2], space...) for t in 0:generations],
  :age => 3,
  :reproduction_start_age => 1,
  :reproduction_end_age => 3,
  :recombination => 1,
  :initial_energy => 0,
  :bottleneck_function => [fill(0.0, space...) for t in 0:generations]
)

## 2. Model parameters
##----------------------------------------------------------------

#NB this dict should be called model_parameters

model_parameters = Dict(
  :species => [species1, species2],
  :generations => generations,  # number of simulation steps
  :space => space,
  :metric => "chebyshev",  # how many neighbors a space site has. "chebyshev" metric means that the r-neighborhood of a position are all positions within the hypercube having side length of 2*floor(r) and being centered in the origin position. "euclidean" metric means that the r-neighborhood of a position are all positions whose cartesian indices have Euclidean distance ≤ r from the cartesian index of the given position.
  :periodic => false,  # whether boundaries of the space are connected
  :resources => [env_resources(x) for x in 0:generations],  # a function that returns a vector of integers for the available resources per site at a given time. It accepts only the model object as the input.
  :interactions => [0.1 0.0; 0.0 -1.0],  # How individuals from different species interact. value  is strength of interaction (between 0 and 1). Sign is the direction of interaction where positive means similar individuals interact more strongly and negative is dissimilar ones tend to interact more. 
  :food_sources => [1.0 0.0; 0.5 0.0],  # What each species feeds on (consumption rate). Has priority over interactions. Non-zero diagonal means the food resource is from the environment. It will be read from rows (species order) to columns (species order).
  :seed => 3
)