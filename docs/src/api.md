# API

## Parameters

The parameters below are required for any simulation. They should be in a dictionary object. The dictionary keys should be the parameters names as `Symbol` (with a colon ":" before the name). See the [Tutorial](@ref) page for an example of the parameters dictionary.

* __m:__ A tuple or an array specifying the ploidy of each species. Currently only support m=1 (haploid) and m=2 (diploid).
* __P:__ A tuple or an array specifying the number of traits _p_ for each species.
* __L:__ A tuple or an array specifying the number of loci _l_ for each species.
* __R:__ A tuple or an array specifying growth for each species. Growth rates are for a logistic growth model, where $N_{t+1} = N_t + r\times N\times (1 - ((N/K))$, where N is population size, t is time, r is growth rate and K is carrying capacity. If _r=0_, population size remains constant.
* __C:__ A matrix containing competition coefficients between each pair of species. A competition coefficient denotes the strength of competition exerted by an individual of species j on an individual of species i. It uses the Lotka-Voltera equation $Ni_{t+1} = Ni_t + r\times N\times (1 - ((Ni + cNj)/K)$ where c is competition coefficient. When competition coefficient is positive, population j competes with population i. If negative, population j helps population i to grow. And if 0, population j does not affect population i. If $c_{ij} > 0$ and $c_{ji} > 0$, both populations are in competition, if $c_{ij} > 0$ and $c_{ji} < 0$, species i is a parasite of species j. If $c_{ij} < 0$ and $c_{ji} < 0$, the two species have a mutualistic relationship. If $c_{ij} < 0$ and $c_{ji} = 0$, they have a commensal relationship.
* __A:__ A tuple or an array of genotype-phenotype matrices, one for each species. Each matrix shows the amount of contribution of each locus to each phenotype. It has P rows and L columns. Make sure that a row is not all zero (a trait is controlled by no locus).
* __Y:__ A tuple or an array  of selection coefficients for each species.
* __T:__ A tuple or an array of arrays, where each inner array specifies optimal phenotypes θ for each species. Each inner array should be of length _p_ (number of traits) of its corresponding species.
* __Ω:__ A tuple or an array of matrices, each of which ω represents a covariance matrix of the selection surface. Each matrix is of size $p\times p$.
* __M:__ A tuple or an array of mutation rates _μ_ for each species.
* __N:__ A dictionary where each key is a node number and its value is a tuple for population size of each species at that node. This dictionary does not have to have a key for all the nodes, but it should have a value for all the species.
* __K:__ A dictionary where each key is a node number and its value is tuple of carrying capacities K of the node for each species. The dictionary should have a key for all the nodes and carrying capacity for each species per node.
* __migration_rates:__ An array of matrices, each matrix shows of migration rates between each pair of nodes for a species. The rows and columns of the matrix are node numbers in order. If instead of a matrix, there is `nothing`, no migration occurs for that species.
* __E:__ A tuple  or an array of the variance of a normal distribution ε representing environmental noise for each species.
* __generations:__ number of generations to run the simulation.
* __space:__ Either a tuple of size 2 or 3 for a grid size or a `SimpleGraph` object for an arbitrary graph. If it is a tuple, a grid is built internally
* __moore:__ Whether nodes in the grid have 8 neighbors (Moore neighborhood). Default is false, i.e. cells only have 4 neighbors.

## Simulation outline

Within each time-step, the following occurs:

1. Mutation
2. Fitness update
3. Migration
4. Reproduction (only for diploids)
5. selection

### Mutation

The genotype vector _y_ and pleiotropy matrix _B_ of each individual mutates.

_y_ mutates by adding random numbers from a normal distribution with mean 0 and standard deviation from the δ to the current values of _y_. δ is specified by the M parameter.

_B_ mutates by randomly switching 0s and 1s with probability given in parameter MB.

### Fitness update

Fitness of each individual updates after mutation. Fitness is $W = exp(γ \times transpose(z - θ)\times inv(ω)\times (z - θ))$, where is the phenotype vector ($z = By + μ$), γ is selection coefficient, θ is optimum phenotypes vector, and ω is covariance matrix of selection surface. 

### Migration

Each agent moves with probabilities given in *migration_rates* to another node.

### Reproduction

When a species is diploid, they sexually reproduce. To that end, individuals of the same species in the same location are randomly paired. Each pair produces one offspring. Then the parents die.

To produce an offspring, each parent contributes to half of the offspring's genotype _y_ and pleiotropy matrix _B_. The genes coming from each parent are randomly assigned.

### Selection

A number of individuals _n_ are selected for the next generation via sampling with replacement weighted by individuals' fitness values. _n_ is calculated using the Lotka-Voltera model for population competition $Ni_{t+1} = Ni_t + r\times N\times (1 - ((Ni + cNj)/K)$ where N is population size, t is time, r is growth rate and K is carrying capacity, and c is competition coefficient. Briefly, each population growth with a logistic model when it is not affected by other species. Otherwise, its growth increases or decreases depending on its interactions with other species.

## Data collection

TODO
