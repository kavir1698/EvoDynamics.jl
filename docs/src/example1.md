```@meta
EditURL = "<unknown>/examples/example1.jl"
```

# Simple Wright-Fisher

In this example, we demonstrate how to create and run simple Wright-Fisher simulations using EvoDynamics.jl. The Wright-Fisher model is a classic population genetics model that describes the genetic drift in a population due to random sampling during reproduction. EvoDynamics.jl allows us to simulate the dynamics of a single haploid species in an unstructured space, with two genes affecting both biotic and abiotic traits.

To get started, we first import the EvoDynamics.jl package:

````@example example1
using EvoDynamics
````
Next, we define the model parameters for our simple one-species Wright-Fisher simulation. The parameters are organized in a Julia (.jl) file as follows:

```julia
generations = 20
space = (1, 1)

## 1. Species parameters

species1 = Dict(
  :name => "a",
  :number_of_genes => 2,
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
  :periodic => false,
  :resources => [fill(200, 1, 1) for i in 0:generations],
  :interactions => [-0.1],
  :food_sources => [1.0],
  :seed => nothing
)
```

In this example, we set the number of generations to 20 and define an unstructured space with a single site (1, 1). The species parameters are specified for a single species named "a", which has two genes influencing two phenotypes (one abiotic and one biotic). The model parameters include information about the species, generation count, space dimensions, metric type, periodic boundaries, resource availability, species interactions, and random seed.

To run the simulation, we use the runmodel function, passing the parameter file path as an argument:

````@example example1
param_file = "../examples/paramfile1.jl"
agentdata, modeldata = runmodel(param_file);
nothing #hide
```

By default, `runmodel` uses the following functions for data collection: `mean_fitness_per_species`, `species_N`. They collect the mean fitness and the population size of each species per time step. This is just a sample data. However, you can customize your own data collection function to gather specific information from the model (see [Collecting data](@ref)).

Finally, we can examine the results of the simulation by inspecting the `modeldata`:

````@example example1
modeldata
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

