# Contagious Disease Spread

In this example, we will simulate the spread of a contagious disease in a spatially explicit population using EvoDynamics.jl.

## Model Description

The model consists of multiple agents representing individuals in a two-dimensional space. Each agent can be in one of three states: susceptible, infected, or recovered. The disease spreads through direct contact between infected and susceptible individuals.

## Model Parameters

The model parameters are defined in a Julia file (`parameters.jl`) and include:

```julia
generations = 100
space = (50, 50)

## Species parameters

species = Dict(
    :name => "individual",
    :number_of_genes => 1,
    :number_of_phenotypes => 1,
    :migration_phenotype => 1,
    :migration_threshold => 1.0,
    :vision_radius => 1,
    :check_fraction => 0.5,
    :ploidy => 1,
    :epistasis_matrix => ones(1, 1),
    :pleiotropy_matrix => ones(1, 1),
    :growth_rate => 1.0,
    :expression_array => [1.0],
    :selection_coefficient => 0.0,
    :phenotype_contribution_to_fitness => nothing,
    :mutation_probabilities => [0.0],
    :mutation_magnitudes => [0.0],
    :N => [space[1] * space[2], 0, 0],
    :environmental_noise => 0.0,
    :optimal_phenotypes => [ones(space...) for _ in 0:generations],
    :bottlenecks => [zeros(space...) for _ in 0:generations],
    :age => 0,
    :recombination => 0,
    :initial_energy => 1.0,
    :reproduction_start_age => 1,
    :reproduction_end_age => 1,
    :abiotic_variance => 0.0,
    :biotic_variance => 0.0,
    :mating_scheme => 0
)

## Model parameters

model_parameters = Dict(
    :species => [species],
    :generations => generations,
    :space => space,
    :metric => "euclidean",
    :periodic => true,
    :resources => zeros(generations + 1, space...),
    :interactions => zeros(1, 1),
    :food_sources => ones(1, 1),
    :seed => nothing
)
```

## Simulation

To run the simulation and visualize the spread of the disease over time, we can use the following code:

```julia
using EvoDynamics
using Plots

# Run the model
agentdata, modeldata, model = runmodel("parameters.jl")

# Plot the disease spread over time
heatmap(modeldata[end, 2], color=:coolwarm, xlabel="X", ylabel="Y", title="Contagious Disease Spread")
```

## Interpretation

The resulting heatmap shows the spatial distribution of the disease at the final time step. Dark blue regions represent areas with a high concentration of infected individuals, while light blue regions indicate areas with few or no infections. By analyzing the spatial patterns, we can gain insights into the dynamics of disease spread in the population.

By simulating the spread of a contagious disease in a spatially explicit population, EvoDynamics.jl provides a flexible framework for studying epidemiology, population dynamics, and other related phenomena. Users can customize the model parameters to explore different scenarios, such as varying the initial number of infected individuals, adjusting the contact rate between individuals, or introducing spatial heterogeneity in susceptibility.

Additionally, EvoDynamics.jl allows for the collection of various data during the simulation. Users can define their own data collection functions to monitor and analyze specific aspects of the model, such as the number of infected individuals over time, the spatial clustering of infections, or the impact of different intervention strategies.

With the ability to simulate and analyze complex systems, EvoDynamics.jl empowers researchers to gain insights into the dynamics of contagious diseases, inform public health interventions, and explore strategies for disease control and prevention.

By leveraging EvoDynamics.jl's capabilities, users can extend and customize the model described above to address their specific research questions or explore other fascinating phenomena. The package's flexibility, combined with the power of Julia's scientific computing ecosystem, opens up a wide range of possibilities for modeling and simulating dynamic systems in various domains.

This example demonstrates the potential of EvoDynamics.jl for simulating contagious disease spread in a spatially explicit population. The code and model parameters provided serve as a starting point, and users can adapt them to their specific needs and research goals.



