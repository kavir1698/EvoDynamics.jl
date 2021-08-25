# API

## Parameters

Parameters of a model should be put in a YAML file with the structure below. Note that spaces and indentations are meaningful in YAML. Indentations should be spaces not tabs.

See [Simple Wright-Fisher](@ref) and [Predator prey](@ref) for complete parameter file examples.

```yml
species:
  1:
    name: a
    parameter 1: ...
    parameter 2: ...
    ...
  2:
    name: b
    parameter 1: ... 
    parameter 2: ... 
    ...
model:
  model parameter 1: ...
  model parameter 2: ...
  ...
```

The file has two main levels: `species` and `model`. `species` stores species specific parameters as many different species as you want.

Since we cannot write a matrix in a YAML file, any parameter that is a matrix should be converted to a vector. In Julia, you can do this by `vec(yourmatrix)`.

### Species specific parameters

Each species should have the following parameters. The order that you write these parameters does not matter.

* __name__: a name for the species.
* __number of genes__: An _integer_ for number of genes that the species has.
* __ploidy__: Either 1 for haploid or 2 for diploid genomes. Diploids can have recombination.
* __number of phenotypes__: An _integer_ for the number of phenotypes that the species has.
* __abiotic phenotypes__: An _array of integers_ (e.g. "[1,2]") specifying abiotic phenotypes among all phenotypes. Abiotic phenotypes determine how the species interacts with the environment.
* __biotic phenotypes__: An _array of integers_ (e.g. "[3]") specifying biotic phenotypes among all phenotypes. Biotic phenotypes determine how the species interacts with other individuals from the same or different species.
* __migration phenotype__: An _integer_ specifying the phenotype that determines migration trait. If the species does not migrate, put 0.
* __migration threshold__: The phenotypic value of the migration phenotype after which an individual can migrate.
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
* __reproduction start age__: The age at which individuals can reproduce.
* __reproduction end age__: The age after which individuals cannot reproduce.
* __recombination__: Mean of a Poisson distributions for number of crossing overs per sexual reproduction.
* __initial energy__: A parameter for parental care of infants. Values more than 0 indicate that newly born individuals can survive for a number of times without requiring food from the environment/other species. The consumption rate (i.e. how many generations this initial energy suffices) is determined by the sum of the corresponding rows in "food sources" model parameter.
* __bottleneck function__: Forcefully kill specific agents. Path to a file that has a function named "bottleneck". The function accepts two arguments: agent and model. It returns true or false for death or survival of the individual, respectively.
* __bottleneck times__: A boolean matrix (zeros and ones) turned into a vector where rows are sites and columns are time steps.  Determines time steps at which the bottleneck function should be activated.

### Model parameters

* __generations__: An _integer_ for the number of steps the model will run.
* __space__: A _vector of two integers_ that determine the size of a grid for space.
* __metric__: Either "chebyshev" or "euclidian". Determines how many neighbors a space site has. "chebyshev" metric means that the r-neighborhood of a position are all positions within the hypercube having side length of 2*floor(r) and being centered in the origin position. "euclidean" metric means that the r-neighborhood of a position are all positions whose cartesian indices have Euclidean distance ≤ r from the cartesian index of the given position.
* __periodic__: _Boolean__ (true or false) to determine whether boundaries of the space are connected or not.
* __resources__: A _vector of integers_ determining available resources (e.g. vegetation) per site per time step.
* __interactions__: A species-species interaction _matrix of numbers_ determining how individuals from different species interact. Each value  is strength of interaction (between 0 and 1). Sign (+/-) is the direction of interaction where positive means similar individuals interact more strongly and negative is dissimilar ones tend to interact more.
* __food sources__: A species-species food _matrix of numbers_ determining what each species feeds on (consumption rate). Non-zero diagonal means the food resource is from the environment. Non-diagonals mean an species (in the rows) feeds on another species (in the columns). Numbers can be zero or any positive number. The magnitude of the number determines how many generations can an individual live off of given one unit of the food source. For example, if a diagonal is 2, it means that the species will eat one unit of the environmental resources and that is enough for it to live two steps.
* __seed__: Either an _integer_ or _Null_ for random number generator seed.

## Simulation outline

The simulations are fully agent-based, meaning that agents do not receive any model-level knowledge for what happens to them. The following steps happen in order to agents that are activated in a random order.

1. Grow one step older.
2. Migrate.
3. Burn energy.
4. Eat from the environment if the species can.
5. Interact with other individuals.
   1. If meeting another individual of the same species but with different sex, try reproduction.
6. If haploid, reproduce.
7. Survive. Agents die if they are too old, or do not have enough energy, or by chance given its fitness.
8. Go through the bottleneck function if it should be activated at the current time.

### Migration

If a species have the capability to migrate (i.e. _migration phenotype_ is not 0), then at each time step, it is checked whether the phenotypic value of the migration phenotype is above _migration threshold_. If it is, then the agent checks a random _check fraction_ of the sites neighboring it and moves to the most suitable site. Suitability is calculated by comparing the agent's abiotic phenotype and the optimal abiotic phenotypic value for each site. If the checked neighbors have worse conditions than the current site, the agent does not move.

### Burn energy

Agent's energy is reduced by one unit.

### Eat

If the agent is able to eat from the environment (the diagonal of _food sources_ at the corresponding row and column is non-zero) and there is any environmental resources left at the agent's site, then its energy level is boosted by as much as the value in the corresponding element at _food sources_.

### Interacting with other individuals

At each time, the agent interacts with as many agents as there are species in the model. If there are fewer agents in the site than number of species, the agent interacts with all of them. Otherwise, it interacts with at most one randomly picked individual from each species. If most of the individuals at the site are from the first species, then it is more likely that the agent interacts with an individual from the first species and with no individual from the second species. This set up allows interactions between species be dependent on the population size of the species at each site.

For each pair of interacting individuals, we first check whether one individual feeds on the other one. If so, the hunt is successful with a probability proportional to the average phenotypic distance between the biotic phenotypes of the two individuals.

Phenotypic distance between two species depends on the sign of the corresponding element in _interactions_ matrix. If it is negative, then the hunt is more successful if the average phenotype of the two individuals is more different.

The phenotypic distance between two individuals is the average distance between all pairs of biotic phenotypes between the two individuals. To calculate the distance between two phenotypic values, we use the following formula:  

$$
| 0.5 - \text{cdf}(\text{Normal}(\text{ph1},  1), \text{ph2})|
$$

, where cdf is cumulative density function, Normal is a normal distribution with mean ph1 (phenotypic value 1) and variance 1, and ph2 is phenotypic value 2.

If the two individuals are not predator-prey, and they are from the same species but different sexes and they are both in reproduction age, then they try to reproduce. Otherwise, if they are from different species and the two species interact with each other (the corresponding values in the _interactions_ matrix is non-zero) they interact with each other. If the individual 1 to 2 has a positive value in the _interactions_ matrix, then it increases the fitness of individual 2.  If the element is zero, it does not increase or decrease the fitness of individual 2. And it the element is negative, it decreases the fitness of individual 2. Similarly, individual 2 can affect the fitness of individual 1. The two interactions do not need to be symmetric. The magnitude of change in the fitness due to interaction is determined by the phenotypic distance between the two individuals and interaction value in the _interactions_ matrix.

### Reproduction

If an individual is haploid and is in reproduction age, it reproduces an identical offspring to itself except that the expression array, pleiotropy matrix, and the epistasis matrix of the offspring will mutate given the probabilities and magnitues in the _mutation probabilities_ and _mutation magnitudes_ matrices. The number of offsprings is a random number from a Poisson distribution with a mean equal to the species' growth rate times the fitness of the individual.

If two diploid individuals are mated to reproduce, their reproduction success is proportional to their phenotypic similarity. Number of offsprings is a random number from a Poisson distribution with mean equal to the reproduction success of the two individuals times the growth rate of the species.

Each offspring of diploid individuals inherits a gamete from each parent. A gamete is a half of the expression array, pleiotropy matrix, and epistasis matrix. If recombination is allowed (> 0), then the gamets undergo crossing over. The number of crossing overs is a random number from a Poission distributioin with mean equal to the _recombination_ parameter.

#### Survival

At each step, the agent may die due to several factors. First, it dies if it has negative energy (has not been able to eat enough). Second, it dies if it is too old (> max age). Third, it dies with a probability negatively correlated to its fitness. This probability is adjusted to the _selection coefficient_. If the selection coefficient is zero, then the individual does not die. If it is one, the fitness determines survival 100%.

#### Fitness

Fitness of individuals is their abiotic fitness, i.e. their abiotic phenotypic distance to the optimal phenotypes in the environment. as such, their fitness can be different at different sites.

#### Mutation

Mutation can happen at three levels: changing the expression of each gene, changing the pleiotropy matrix, and changing the epistasis interactions between genes. The probability that a mutation occurs at each of these levels is controlled by parameter _mutation probabilities_. And size of mutations when they occur are controlled by the _mutation magnitudes_ parameters, which are the variances of a normal distribution. A mutation in gene expression and epistasis matrix is carried by adding a random number from a normal distribution to the existing values. The pleiotropy matrix is mutated by switching 0s and 1s with the probability given in the third element of _mutation magnitudes_.

## Data collection

The interface to the model is from the `runmodel` function (see [Tutorial](@ref)).

EvoDynamics.jl uses [Agents.jl](https://github.com/JuliaDynamics/Agents.jl) underneath. See [Agents.jl&#39;s documentation](https://juliadynamics.github.io/Agents.jl/dev/) for details about writing functions to collect any data during simulations. Here, we explain the specific implementation of the model.

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
