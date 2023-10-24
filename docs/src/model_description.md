# Model Parameters and Simulation Outline

## Parameters

All model parameters are stored in a Julia (.jl) file and divided into two sections: species-specific parameters and model-specific parameters. Parameters for each species and the model are written in separate Dict objects, with keys represented as Symbols (using a colon : before the names).

To parameterize a model, you can copy an existing set of parameters from the Examples section and modify them accordingly.

### 1. Species specific parameters

Each species should have a dictionary object containing the following parameters. The order in which they are written does not matter.

* __name__: A `String` representing the name of the species.
* __number\_of\_genes__: An _integer_ indicating the number of genes that the species has.
* __ploidy__: Either 1 for haploid or 2 for diploid genomes. Diploids may undergo recombination.
* __number\_of\_phenotypes__: An _integer_ indicating the number of phenotypes that the species has.
* __abiotic\_phenotypes__: An _array of integers_ specifying abiotic phenotypes among all phenotypes. Abiotic phenotypes determine how the species interacts with the environment.
* __biotic\_phenotypes__: An _array of integers_ specifying biotic phenotypes among all phenotypes. Biotic phenotypes determine how the species interacts with other individuals from the same or different species.
* __migration\_phenotype__: An _integer_ specifying the phenotype that corresponds to the migration trait. If the species does not migrate, put 0.
* __migration\_threshold__: The phenotypic value of the migration phenotype after which an individual can migrate.
* __vision\_radius__: A _number_ determining the radius of neighboring sites that the agent can see before migration.
* __check\_fraction__: A _floating-point number_ between 0 and 1 representing the fraction of visible sites that the agent can check and decide whether to migrate to.
* __epistasis\_matrix__: An epistasis matrix of size `l x l`, where `l` is the product of the number of genes and ploidy. The epistasis matrix specifies the direction (positive or negative) and size of the effect of one locus on other loci. For example, if at row 1 and column 2 is a value 0.2, it means that locus 1 affects locus 2 by increasing the effect of locus 2 (because its positive) with 20% of the effect of locus 1.
* __pleiotropy\_matrix__: A _binary matrix_ (0s and 1s) with size _number of phenotypes_ x _l_. The pleiotropy matrix specifies the phenotypes that each locus affects.
* __expression\_array__: A _vector_ of size _l_ representing the expression amount of each locus determining its effect size.
* __growth\_rate__: The mean of a Poisson distribution for the number of offspring per reproduction. This number is fully realized when the fitness of a haploid individual is 1 or the distance between the biotic phenotypes of two diploid individuals is 0. Otherwise, the actual growth rate is a fraction of this value.
* __selection\_coefficient__: A number between 0 and 1 that determines the importance of fitness. 0 represents a model without selection. To specify time-variable selection coefficient, provide a vector with a length equal to the number of generations plus 1 (to account for generation zero).
* __phenotype\_contribution\_to\_fitness__: The relative contribution of each __abiotic phenotype__ to fitness. This parameter can be `nothing` to denote equal contribution of all abiotic phenotypes, or it can be a vector of numbers with as many elements as there are abiotic phenotypes.
* __mutation\_probabilities__: A _vector of three numbers_ each of which specifies the probability for a different type of mutations: mutation probability of the `expression_array`, `pleiotropy_matrix`, and `epistasis_matrix`, respectively.
* __mutation\_magnitudes__: A _vector of numbers_ with the same size as `mutation_probabilities` that determines the magnitude of mutation for each of the three categories. Specifically, the numbers represent the variances of normal distributions with mean 0 for the expression array and epistasis matrices, and the probability of changing a 0 and 1 in the pleiotropy matrix. When mutating an offspring, it first checks whether there will be a mutation (`mutation_probabilities`), and if positive, a mutation with a magnitude determined by _mutation\_magnitudes_ occurs.
* __N__: A _vector of integers_ representing the initial number of individuals at each site.
* __environmental\_noise__: A number for the variance of a normal distribution with mean 0 that will be added to the phenotypes.
* __optimal\_phenotypes__: A _vector of matrices_ representing the optimal phenotypes at each time step. Each matrix represents the optimal phenotypes for each site in space. The optimal phenotype at each site is a vector as long as the `abiotic_phenotypes`. The vector is as long as the number of model steps (generations) plus 1 (for time zero).
* __age__: An _integer_ for the maximum age of individuals of this species.
* __reproduction\_start\_age__: The age at which individuals can reproduce.
* __reproduction\_end\_age__: The age after which individuals cannot reproduce.
* __recombination__: The mean of a Poisson distribution for the number of crossing overs per sexual reproduction.
* __initial\_energy__: A parameter for parental care of infants. Values greater than 0 indicate that newly born individuals can survive for a certain number of time steps without requiring food from the environment or other species. The consumption rate is always constant for all species at one unit per time step. Use with caution as having an initial energy larger than zero can lead to infinite population growth when individuals without food reproduce and their offspring also reproduce without food. For example, this happens when start age of reproduction is 1 and initial energy is larger than 0.
* __bottlenecks__: A vector of matrices, each of which has the same size as the space and contains the probability of death at each site. The vector is as long as the number of model steps (generations) plus 1 (for time zero). Use this if you want to impose a killing event in the model at specific times and places. Otherwise, use 0.0 for all the values in the matrices.
* __mating\_scheme__: One of the following: -1, 0, 1. Only used in sexual reproduction. Zero means randomly paired couples are likely to have as many offspring as any other pair. Their phenotypes do not affect their reproductive success. -1 represents disassortative mating, where the more different the phenotypes of the pair, the more children they have. 1 represents assortative mating, where individuals with similar phenotypes are more likely to reproduce.
* __abiotic\_variance__: The variance of a normal distribution used in determining the phenotypic distance of agents to the optimal environmental phenotypes. The larger the variance, the less important is the distance.
* __biotic\_variance__: The variance of a normal distribution used in determining the biotic phenotypic distance between two individuals (used in any kind of interaction). Larger values mean that all pairs are equally likely to interact, regardless of their phenotype difference.

### 2. Model parameters

This dictionary should be named "model_parameters".

* __species__: A vector containing the dictionaries of the species.
* __generations__: An integer indicating the number of steps the model will run.
* __space__:A tuple of two integers (e.g., `(3,2)`) determining the size of space. If you do not want any spatial structure, use `(1,1)`.
* __metric__: Either "chebyshev" or "euclidean". Determines how many neighbors a space site has. "chebyshev" metric means that the r-neighborhood of a position includes all positions within the hypercube with a side length of 2 * floor(r) and centered at the origin position. "euclidean" metric means that the r-neighborhood of a position includes all positions whose Cartesian indices have Euclidean distance ≤ r from the Cartesian index of the given position.
* __periodic__: A Boolean value (true or false) determining whether the boundaries of the space are connected or not.
* __resources__: A vector of matrices where each matrix has the size of the space and determines the amount of resources per site. The vector is as long as the number of model steps (generations) plus 1 (for time zero). The elements should be floating point numbers.
* __interactions__: A species-species interaction matrix of floating-point numbers determining how individuals from different species interact. Each value represents the strength of interaction (between 0 and 1). The sign (+/-) indicates the direction of interaction, where a positive value means similar individuals interact more strongly, and a negative value means dissimilar individuals tend to interact more.
* __food\_sources__: A species-species food matrix of floating-point numbers determining what each species feeds on (consumption rate). Non-zero values on the diagonal indicate that the food resource is from the environment. Off-diagonal values indicate that a species (in the rows) feeds on another species (in the columns). Numbers can be zero or any positive number. The magnitude of the number determines how much energy the predator gains by eating the prey. All species use 1 unit of energy per time step.
* __seed__: Either an integer or `nothing` for a random number.


## Simulation outline

The simulations are fully agent-based, meaning that agents do not receive any model-level knowledge about what happens to them. The following steps occur in order for agents that are activated in a random order:

1. Grow one step older.
2. Migrate.
3. Burn energy.
4. Eat from the environment if the species can.
5. Interact with other individuals.
6. Reproduce.
7. Survive. Agents die if they are too old, or do not have enough energy, or by chance given its fitness.
8. Kill agents with a probability given in `bottlenecks`.

### Migration

If a species has the capability to migrate (i.e., `migration_phenotype` is not 0), then at each time step, it is checked whether the phenotypic value of the migration phenotype is above the `migration_threshold`. If it is, then the agent checks a random `check_fraction` of the neighboring sites and moves to the most suitable site. Suitability is calculated by comparing the agent's abiotic phenotype and the optimal abiotic phenotypic value for each site. If the checked neighbors have worse conditions than the current site, the agent does not move.


### Burn energy

Agent's energy is reduced by one unit.

### Eat

If the agent is able to eat from the environment (the diagonal of `food_sources` at the corresponding row and column is non-zero) and there are any environmental resources left at the agent's site, its energy level is boosted by the value in the corresponding element at `food_sources`.


### Interacting with other individuals

At each time step, the agent interacts with individuals from each species in the model. If there are fewer agents at the site than the number of species, the agent interacts with all of them. Otherwise, it interacts with a random selection of individuals in the regions. If most of the individuals at the site are from the first species, then the agent is more likely to interact with an individual from the first species and with no individual from the second species. This setup allows interactions between species to depend on the population size of the species at each site. At each time step, an individual from one species can interact with individuals from each other species at most once.

For each pair of interacting individuals, it is first checked whether one individual feeds on the other. If so, the hunt is successful with a probability proportional to the average phenotypic distance between the biotic phenotypes of the two individuals.

The phenotypic distance between two species depends on the sign of the corresponding element in the `interactions` matrix. If it is negative, the hunt is more successful when the average phenotype of the two individuals is more different.

The phenotypic distance between two individuals is calculated as the average distance between all pairs of biotic phenotypes between the two individuals. The distance between two phenotypic values is determined using the formula:

$ | 0.5 - \text{cdf}(\text{Normal}(\text{ph1},  1), \text{ph2})| $

where `cdf` is the cumulative density function, `Normal` is a normal distribution with mean `ph1` (phenotypic value 1) and variance 1, and `ph2` is phenotypic value 2.

If the two individuals are not predator-prey and they are from the same species but different sexes, and both are in reproductive age, they are marked to reproduce in the next stage. Otherwise, if they are from different species and the two species interact with each other (the corresponding values in the `interactions` matrix are non-zero), they interact with each other. If the individual from species 1 to species 2 has a positive value in the `interactions` matrix, it increases the fitness of the individual from species 2. If the element is zero, it does not affect the fitness of the individual from species 2. And if the element is negative, it decreases the fitness of the individual from species 2. Similarly, the individual from species 2 can affect the fitness of the individual from species 1. The two interactions do not need to be symmetric. The magnitude of change in fitness due to interaction is determined by the phenotypic distance between the two individuals and the interaction value in the `interactions` matrix.

### Reproduction

If an individual is haploid and in reproductive age, it reproduces an identical offspring to itself, except that the expression array, pleiotropy matrix, and the epistasis matrix of the offspring will mutate based on the probabilities and magnitudes in the `mutation_probabilities` and `mutation_magnitudes` matrices. The number of offspring is a random number from a Poisson distribution with a mean equal to the species' growth rate times the fitness of the individual.

If two diploid individuals mate to reproduce, their reproductive success is proportional to their phenotypic distance (depending on `mating_scheme`). The number of offspring is a random number from a Poisson distribution with a mean equal to the reproductive success of the two individuals times the growth rate of the species.

Each offspring of diploid individuals inherits a gamete from each parent. A gamete is half of the expression array, pleiotropy matrix, and epistasis matrix. If recombination is allowed (value > 0), then the gametes undergo crossing over. The number of crossing overs is a random number from a Poisson distribution with a mean equal to the `recombination` parameter.

#### Survival

At each step, the agent may die due to several factors. First, it dies if it has negative energy (has not been able to eat enough). Second, it dies if it is too old (greater than the maximum age). Third, it dies with a probability negatively correlated to its fitness. This probability is adjusted by the `selection_coefficient`. If the selection coefficient is zero, then the individual does not die. If it is one, the fitness determines a 100% survival rate.

#### Fitness

The fitness of individuals is determined by their abiotic fitness, i.e., their abiotic phenotypic distance to the optimal phenotypes in the environment. Thus, their fitness can vary at different sites.

#### Mutation

Mutation can occur at three levels: changing the expression of each gene, changing the pleiotropy matrix, and changing the epistasis interactions between genes. The probability that a mutation occurs at each of these levels is controlled by the `mutation_probabilities` parameter. The size of mutations, when they occur, is controlled by the `mutation_magnitudes` parameter, which represents the variances of normal distributions. A mutation in the gene expression and epistasis matrix is carried out by adding a random number from a normal distribution to the existing values. The pleiotropy matrix is mutated by switching 0s and 1s with the probability given in the third element of `mutation_magnitudes`.
