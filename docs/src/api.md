# API

## Parameters

The parameters below are required for any simulation. They should be in a dictionary object. The dictionary keys should be the parameters names as `Symbol` (with a colon ":" before the name). See the [Tutorial](@ref) page for an example of the parameters dictionary.

* __ploidy:__ A tuple specifying the ploidy of each species. Currently only support ploidy=1 (haploid) and ploidy=2 (diploid).
* __phenotypes:__ A tuple specifying the number of phenotypes _p_ for each species.
* __ngenes:__ A tuple specifying the number of genes _l_ for each species.
* __growthrates:__ A tuple specifying growth rate for each species. Growth rates are for a logistic growth model, where $N_{t+1} = N_t + r\times N_t \times (1 - ((N_t / K))$, where N is population size, t is time, r is growth rate and K is carrying capacity. If _r=0_, population size remains constant.
* __interactionCoeffs:__ A matrix containing competition coefficients between each pair of species. A competition coefficient $c_{ij}$ denotes the strength of competition exerted by an individual of species j on an individual of species i.
* __interaction_equation__: either "lotkaVoltera\_generalized" (if you want host-parasite or predator prey dynamics) or "lotkaVoltera\_competition". Default is "lotkaVoltera\_competition". The LV competition equation is as follows. It uses the Lotka-Voltera equation $Ni_{t+1} = Ni_t + r\times Ni_t \times (1 - ((Ni_t + \sum_{j=1}^{Nspecies} c_{ij}Nj_t)/Ki)$ where c is competition coefficient. When competition coefficient is positive, population j competes with population i. If negative, population j helps population i to grow. And if 0, population j does not affect population i. If $c_{ij} > 0$ and $c_{ji} > 0$, both populations are in competition, if $c_{ij} > 0$ and $c_{ji} < 0$, species i is a parasite of species j. If $c_{ij} < 0$ and $c_{ji} < 0$, the two species have a mutualistic relationship. If $c_{ij} < 0$ and $c_{ji} = 0$, they have a commensal relationship. It can also be `nothing`. The generalized Lotka-Voltera equation is $ x_{t+1} = x_t + D(x_t)(r + Cx) $, where $ x $ is an array population densities of all species ($ N/K $), $ r $ is an array of growth rates when growing alone, C is interaction coefficient matrix.
* __pleiotropyMat:__ A tuple of pleiotropy matrices, one for each species. A pleitropy matrix is a binary matrix with size $p \times l$ that specifies the phenotypes that each locus affects.
* __epistasisMat:__ A tuple of epistasis matrices, one for each species. An epistasis matrix is of size $l \times l$ and specifies the direction (positive or negative) and size of effect of one loci on the other loci. For example, if at row 1 and column 2 is a value 0.2, it means that locus 1 affects locus 2 by increasing the effect of locus 2 (because its positive) with 20% of the effect of locus 1. The effect of loci are themselves specified in the $q$ vector.
* __expressionArrays:__ A tuple of expression arrays, one for each species. An expression array $q$ shows the effect size of a locus. It can be thought of expression level of a gene.
* __selectionCoeffs:__ A tuple  of selection coefficients for each species.
* __optPhenotypes:__ A tuple of arrays, where each inner array specifies optimal phenotypes θ for each species. Each inner array should be of length _p_ (number of phenotypes) of its corresponding species.
* __covMat:__ A tuple of matrices, each of which ω represents a covariance matrix of the selection surface. Each matrix is of size $p\times p$.
* __mutProbs:__ A tuple of arrays, one per species. Each inner array specifies the probabilities for different type of mutations. An inner array has three values with the following order: first, the probability that an individual's gene expressions (array $q$) mutate per generation. Second, the probability that an individual's pleiotropy matrix mutates per generation. Third, the probability that an individual's epistasis matrix mutates per generation. If an individual is going to receive a mutation at any of the mentioned levels, the amount of change that it receives is determined by $mutMagnitudes$.
* __mutMagnitudes:__ A tuple of arrays, one for each species. Each inner array has three numbers specifying the amount of change per mutation in gene expression ($q$), pleiotropy matrix ($b$), and epistasis matrix $a$, respectively. The first number is the variance of a normal distribution with mean zero. If an individual's gene expression is going to mutate, random numbers from that distribution are added to the expression of each locus. The second number is a probability that any one element in the pleiotropy matrix will switch - if on, becomes off, and vice versa. The third number is again the variance of another normal distribution with mean zero. Random numbers taken from such distribution are added to each element of the epistasis matrix.
* __N:__ A dictionary where each key is a node number and its value is a tuple for population size of each species at that node. This dictionary does not have to have a key for all the nodes, but it should have a value for all the species.
* __K:__ A dictionary where each key is a node number and its value is tuple of carrying capacities K of the node for each species. The dictionary should have a key for all the nodes and carrying capacity for each species per node.
* __migration_rates:__ An array of matrices, each matrix shows migration rates between each pair of nodes for a species. The rows and columns of the matrix are node numbers in order. Values are read column-wise: each column shows out-migration to all other nodes from a node. If instead of a matrix, there is `nothing`, no migration occurs for that species.
* __E:__ A tuple  of the variance of a normal distribution ε representing environmental noise for each species.
* __generations:__ number of generations to run the simulation.
* __space:__ (default=`nothing`) Either a tuple of size 2 or 3 for a grid size or a `SimpleGraph` object for an arbitrary graph. If it is a tuple, a grid is built internally
* __moore:__ (default=`false`) Whether nodes in the grid have 8 neighbors (Moore neighborhood). Default is false, i.e. cells only have 4 neighbors.
* __periodic:__ (default=`false`) If `space` is 2D, should the edges connect to the opposite side?
* __seed:__ (default=`0`). Seed for random number generator. Only set if >0.

## Simulation outline

Within each time-step, the following occurs:

1. Mutation
2. Fitness update
3. Migration
4. Reproduction (only for diploids)
5. selection

### Mutation

Mutation can happen at three levels: changing the expression of each gene expressionArrays, changing the pleiotropy matrix pleiotropyMat, and changing the epistatic interactions between genes. The probability that a mutation occurs at each of these levels is controlled by parameter mutProbs. And size of mutations when they occur are controlled by parameter mutMagnitudes.
The genotype vector _y_ and pleiotropy matrix _B_ of each individual mutates.

Epistatic matrix _A_ and expression vectors _Q_ mutate by adding their values to random numbers from a normal distribution with mean 0 and standard deviation given in parameter mutMagnitudes.

_B_ mutates by randomly switching 0s and 1s with probability given in parameter mutMagnitudes.

### Fitness update

Fitness of each individual updates after mutation. Fitness is $W = exp(γ \times transpose(z - θ)\times inv(ω)\times (z - θ))$, where is the phenotype vector ($z = pleiotropyMat(Aq) + μ$), γ is selection coefficient, θ is optimum phenotypes vector, and ω is covariance matrix of selection surface. 

### Migration

Each agent moves with probabilities given in *migration_rates* to other nodes.

### Reproduction

When a species is diploid, they sexually reproduce. To that end, individuals of the same species in the same location are randomly paired. Each pair produces one offspring. Then the parents die.

To produce an offspring, each parent contributes to half of the offspring's genotype _y_ and pleiotropy matrix _B_. The genes coming from each parent are randomly assigned.

### Selection

A number of individuals _n_ are selected for the next generation via sampling with replacement weighted by individuals' fitness values. _n_ is calculated using the Lotka-Voltera model for population competition $Ni_{t+1} = Ni_t + r\times N\times (1 - ((Ni + cNj)/K)$ where N is population size, t is time, r is growth rate and K is carrying capacity, and c is competition coefficient. Briefly, each population growth with a logistic model when it is not affected by other species. Otherwise, its growth increases or decreases depending on its interactions with other species.

## Data collection

The interface to the model is from the `runmodel` function (see [Tutorial](@ref)).


EvoDynamics.jl uses [Agents.jl](https://github.com/JuliaDynamics/Agents.jl) underneath. See [Agents.jl's documentation](https://juliadynamics.github.io/Agents.jl/dev/) for details about writing functions to collect any data during simulations. Here, we explain the specific implementation of the model.

There are two main objects from which you can collect data: and agent object of type `AbstractAgent` and a model object of type `ABM`. Both of these types are defined the `Agents.jl` package.

Agent object has the following fields: `id`, `positions`, `species`, `epistasisMat` (epistasis matrix), `pleiotropyMat` (pleiotropy matrix), and `q` (gene expression array).

The model object has the following fields: `space` which is a `GraphSpace` or `GridSpace` object from `Agents.jl`, `agents` that is an array holding all agents, and `properties` which is a dictionary holding all the parameters passed to the model.

To collect data, provide a dictionary where the keys are either agent fields, or `:model`. The value of a key is an array of any number of functions.

If a key is an agent field, all the value of the field from all agents are collected and then aggregated with the functions in the value. For example, to collect mean and median fitness of individuals which is in field `W`, your dictionary will be Dict(:W => [mean, median]).

If a key is `:model`, functions in its value array should be functions that accept a single argument, the model object, and return a single number or a tuple of numbers. For example, this is the default dictionary and its function:

```jl
collect = Dict(:model => [mean_fitness_per_species])

"Returns a tuple whose entries are the mean fitness of each species."
function mean_fitness_per_species(model::ABM)
  nspecies = length(model.properties[:nphenotypes])
  mean_fitness = Array{Float32}(undef, nspecies)
  for species in 1:nspecies
    fitness = mean([i.W for i in values(model.agents) if i.species == species])
    mean_fitness[species] = fitness
  end

  return Tuple(mean_fitness)
end
```
