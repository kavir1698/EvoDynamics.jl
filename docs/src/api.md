# API

## Parameters

Parameters of a model should be put in a YAML file with the structure below. Note that spaces and indentations are meaningful in YAML. Indentations should be spaces not tabs.

See [Simple Wright-Fisher](@ref) and [Predator prey](@ref) for complete parameter file examples.

```yml
- species:
  - id: 1
    parameter 1: ...
    parameter 2: ...
    ...
  - id: 2
    parameter 1: ... 
    parameter 2: ... 
    ...
- model:
  model parameter 1: ...
  model parameter 2: ...
  ...
```

The file has two main levels: `species` and `model`. `species` stores species specific parameters as many different species as you want, starting with `id` 1. 

Since we cannot write a matrix in a YAML file, any parameter that is a matrix should be converted to a vector. In Julia, you can do this by `vec(yourmatrix)`.

### Species specific parameters

Each species should have the following parameters. The order that you write these parameters does not matter.

* __number of genes__: An _integer_ for number of genes that the species has.
* __ploidy__: Either 1 for haploid or 2 for diploid genomes.
* __number of phenotypes__: An _integer_ for the number of phenotypes that the species has.
* __abiotic phenotypes__: An _array of integers_ (e.g. "[1,2]") specifying abiotic phenotypes among all phenotypes. Abiotic phenotypes determine how the species interacts with the environment.
* __biotic phenotypes__: An _array of integers_ (e.g. "[3]") specifying biotic phenotypes among all phenotypes. Biotic phenotypes determine how the species interacts with other individuals from the same or different species.
* __migration phenotype__: An _integer_ specifying the phenotype that determines migration trait. If the species does not migrate, put 0.
* __vision radius__: A _number_ determining the radius of neighboring sites that the agent can see before migration.
* __check fraction__: A _number_ between 0 and 1 showing the fraction of the visible sites to the agent that it can check and decide whether to migrate to.
* __epistasis matrix__: An epistasis matrix is of size $l \times l$, where _l_ is the product of _number of genes_ and _ploidy_. Epistasis matrix specifies the direction (positive or negative) and size of effect of one locus on other loci. For example, if at row 1 and column 2 is a value 0.2, it means that locus 1 affects locus 2 by increasing the effect of locus 2 (because its positive) with 20% of the effect of locus 1.
* __pleiotropy matrix__: A _binary matrix_ (0s and 1s) with size _number of phenotypes_ times _l_. The pleiotropy matrix specifies the phenotypes that each locus affects.
* __expression array__: A _vector_ of size _l_ that represent the expression amount of each locus determining its effect size.
* __growth rate__: Mean of a Poisson distribution for number of offsprings per reproduction. This number is the maximum mean when fitness of a haploid individual is 1, or the distance between the biotic phenotypes of two diploid individuals is 0.
* __selection coefficient__: A number between 0 and 1 that determines the importance of fitness. 0 would be a model without selection.
* __mutation probabilities__: A _vector of three numbers_ each of which specifies the probability for a different type of mutations: mutation probability of the _expression array_, _pleiotropy matrix_, and _epistasis matrix_, respectively. 
* __mutation magnitudes__: A _vector of numbers_ with the same size as _mutation probabilities_ that determines the magnitude of mutation for each of the three categories. Specifically, the numbers are the variances of normal distributions with mean 0 for expression array and epistasis matrices, and probability of changing a 0 and 1 in in the pleiotropy matrix.
* __N__: A _vector of integers_ for the initial number of individuals at each site.
* __environmental noise__: A number for the variance of a normal distribution with mean 0 that will be added to the phenotypes.
* __optimal phenotype values__: One or more vectors of numbers. Each vector should start with a dash and one indentation level. Each vector is the optimal phenotypes for each site for all abiotic traits. There are as many element as number of sites times number of abiotic traits. The first N elements are for the first abiotic trait, where N is the number of sites, and so on. 
* __optimal phenotypes__: A _vector of integers_ for optimal phenotype indices for each generation, including generation zero.
* __age__: An _integer_ for maximum age of individuals of this species.
* __recombination__: Mean of a Poisson distributions for number of crossing overs per sexual reproduction.
* __initial energy__: A parameter for parental care of infants. Values more than 0 indicate that newly born individuals can survive for a number of times without requiring food from the environment/other species. The consumption rate (i.e. how many generations this initial energy suffices) is determined by the sum of the corresponding rows in "food sources" model parameter.

### Model parameters

* __generations__: An _integer_ for the number of steps the model will run.
* __space__: A _vector of two integers_ that determine the size of a grid for space.
* __metric__: Either "chebyshev" or "euclidian". Determines how many neighbors a space site has. "chebyshev" metric means that the r-neighborhood of a position are all positions within the hypercube having side length of 2*floor(r) and being centered in the origin position. "euclidean" metric means that the r-neighborhood of a position are all positions whose cartesian indices have Euclidean distance ≤ r from the cartesian index of the given position.
* __periodic__: _Boolean__ (true or false) to determine whether boundaries of the space are connected or not.
* __resources__: A _vector of integers_ determining available resources (e.g. vegetation) per site per time step. 
* __interactions__: A species-species interaction _matrix of numbers_ determining how individuals from different species interact. Each value  is strength of interaction (between 0 and 1). Sign (+/-) is the direction of interaction where positive means similar individuals interact more strongly and negative is dissimilar ones tend to interact more.
* __food sources__: A species-species food _matrix of numbers_ determining what each species feeds on (consumption rate). Has priority over interactions. Non-zero diagonal means the food resource is from the environment. It will be read from rows (species order) to columns (species order). 
* __seed__: Either an _integer_ or _Null_ for random number generator seed.

## Simulation outline

The simulations are fully agent-based, meaning that agents do not receive any model-level knowledge for what happens to them. The following steps happen in order to agents that are activated in a random order.

1. Grow one step older.
2. Migrate.
3. Eat basic food if the species can.
4. Consume energy.
5. Interact with other species.
   1. If meeting another individual of the same species but with different sex, reproduce.
6. If haploid, reproduce.
7. Survive if not tool old, has energy and by chance given its fitness.

### Mutation

Mutation can happen at three levels: changing the expression of each gene expressionArrays, changing the pleiotropy matrix pleiotropyMat, and changing the epistasis interactions between genes. The probability that a mutation occurs at each of these levels is controlled by parameter _mutation probabilities_. And size of mutations when they occur are controlled by parameter mutation magnitudes.
The genotype vector _y_ and pleiotropy matrix _B_ of each individual mutates.

Epistasis matrix _A_ and expression vectors _Q_ mutate by adding their values to random numbers from a normal distribution with mean 0 and standard deviation given in parameter mutMagnitudes.

_B_ mutates by randomly switching 0s and 1s with probability given in parameter mutMagnitudes.

### Fitness update

TODO

### Migration

TODO

### Reproduction

TODO

### Survival

TODO

## Data collection

The interface to the model is from the `runmodel` function (see [Tutorial](@ref)).


EvoDynamics.jl uses [Agents.jl](https://github.com/JuliaDynamics/Agents.jl) underneath. See [Agents.jl's documentation](https://juliadynamics.github.io/Agents.jl/dev/) for details about writing functions to collect any data during simulations. Here, we explain the specific implementation of the model.

There are two main objects from which you can collect data: and agent object of type `AbstractAgent` and a model object of type `ABM`. Both of these types are defined the `Agents.jl` package.

Agent object has the following fields that are individual specific: `id`, `pos`, `species`, `epistasisMat` (epistasis matrix), `pleiotropyMat` (pleiotropy matrix), `q` (gene expression array), `biotic_phenotype`, `abiotic_phenotype`, `age`, `sex`, `energy`, `interaction_history`, and `W` (fitness).

The model object has the following fields that can be accessed with the `.` syntax and describe properties of species or the model: `ngenes`, `nphenotypes`, `growthrates`, `selectionCoeffs`, `ploidy`, `optvals` (optimal values), `optinds` (optval indices per generation), `mutProbs` (mutation probabilities), `mutMagnitudes` (mutation magnitudes), `N`, `E` (environmental noise), `generations`, `nspecies`, `migration_traits`, `vision_radius`, `check_fraction`, `migration_thresholds`, `step`, `biotic_phenotypes` (indices of biotic phenotypes per species), `abiotic_phenotypes`, `max_ages`, `food_sources`, `interactions`, `resources`, `recombination`, `initial_energy`.

You can collect data from agents and/or from the model object. To collect data from agents, use the `adata` keyword argument in the `runmodel` function, and to collect data from the model, use the `mdata` keyword. A complete description of the values these keywords take are at [data collection section of the Agents.jl package](https://juliadynamics.github.io/Agents.jl/stable/tutorial/#.-Collecting-data).

For example, we use the function below to count the number of individual per species:

```jl
"Returns the size of each species."
function species_N(model::ABM)
  allagents = model.agents
  if length(allagents) == 0
    return fill(0, model.nspecies)
  else
    counts = countmap([a.species for a in values(model.agents)])
    output = fill(0, model.nspecies)
    for (k, v) in counts
      output[k] = v
    end
    return output
  end
end

using EvoDynamics

agentdata, modeldata, model = runmodel("parameters.yml", mdata=[species_N])
```
