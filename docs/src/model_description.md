# Model description

## Parameters

All parameters of a model are written in a julia (`.jl`) file. Parameters are divided into two sections: species specific and model specific. Parameters for each species and for the model are written in separate dictionary (`Dict`) objects. The keys of the dictionary should be Symbols (use colon : before the names).

An easy way to parameterize a model is to copy an existing set of parameters from the Examples section and modify them.

### 1. Species specific parameters

Each species should have the following parameters in a dictionary objection. The order that you write these parameters does not matter.

* __name__: A `String` as the name of the species.
* __number\_of\_genes__: An _integer_ for number of genes that the species has.
* __ploidy__: Either 1 for haploid or 2 for diploid genomes. Diploids may undergo recombination.
* __number\_of\_phenotypes__: An _integer_ for the number of phenotypes that the species has.
* __abiotic\_phenotypes__: An _array of integers_ (e.g. "[1,2]") specifying abiotic phenotypes among all phenotypes. Abiotic phenotypes determine how the species interacts with the environment.
* __biotic\_phenotypes__: An _array of integers_ (e.g. "[3]") specifying biotic phenotypes among all phenotypes. Biotic phenotypes determine how the species interacts with other individuals from the same or different species.
* __migration\_phenotype__: An _integer_ specifying the phenotype that corresponds to the migration trait. If the species does not migrate, put 0.
* __migration\_threshold__: The phenotypic value of the migration phenotype after which an individual can migrate.
* __vision\_radius__: A _number_ determining the radius of neighboring sites that the agent can see before migration.
* __check\_fraction__: A _floating number_ between 0 and 1 showing the fraction of the visible sites to the agent that it can check and decide whether to migrate to.
* __epistasis\_matrix__: An epistasis matrix is of size $l \times l$, where _l_ is the product of _number of genes_ and _ploidy_. Epistasis matrix specifies the direction (positive or negative) and size of effect of one locus on other loci. For example, if at row 1 and column 2 is a value 0.2, it means that locus 1 affects locus 2 by increasing the effect of locus 2 (because its positive) with 20% of the effect of locus 1.
* __pleiotropy\_matrix__: A _binary matrix_ (0s and 1s) with size _number of phenotypes_ times _l_. The pleiotropy matrix specifies the phenotypes that each locus affects.
* __expression\_array__: A _vector_ of size _l_ that represent the expression amount of each locus determining its effect size.
* __growth\_rate__: Mean of a Poisson distribution for number of offsprings per reproduction. This number is fully realized when fitness of a haploid individual is 1, or the distance between the biotic phenotypes of two diploid individuals is 0. Otherwise, actual growth rate is a fraction of this value. 
* __selection\_coefficient__: A number between 0 and 1 that determines the importance of fitness. 0 is be a model without selection. If you want time-variable selection coefficient, provide a vector as long as number of generations plus 1 (to account for generation zero).
* __phenotype\_contribution\_to\_fitness__: The relative contribution of each __abiotic phenotype__ to fitness. This parameter can be `nothing` to denote equal contribution of all abiotic phenotypes, or it can be a vector of numbers with as many elements as there are abiotic phenotypes.
* __mutation\_probabilities__: A _vector of three numbers_ each of which specifies the probability for a different type of mutations: mutation probability of the _expression\_array_, _pleiotropy\_matrix_, and _epistasis\_matrix_, respectively.
* __mutation\_magnitudes__: A _vector of numbers_ with the same size as _mutation\_probabilities_ that determines the magnitude of mutation for each of the three categories. Specifically, the numbers are the variances of normal distributions with mean 0 for expression array and epistasis matrices, and probability of changing a 0 and 1 in in the pleiotropy matrix. When mutating an offspring, it first checks whether there will be a mutation (_mutation\_probabilities_), and if positive, a mutation with a magnitude determined by _mutation\_magnitudes_ occurs.
* __N__: A _vector of integers_ for the initial number of individuals at each site.
* __environmental\_noise__: A number for the variance of a normal distribution with mean 0 that will be added to the phenotypes.
* __optimal\_phenotypes__: A _vector of matrices_. The vector holds optimal phenotypes at each time step. Each matrix, is the optimal phenotypes for each site in space. The optimal phenotype at each site is a vector as long as `abiotic_phenotypes`. The vector is as long as model steps (generations) plus 1 (for time zero).
* __age__: An _integer_ for maximum age of individuals of this species.
* __reproduction\_start\_age__: The age at which individuals can reproduce.
* __reproduction\_end\_age__: The age after which individuals cannot reproduce.
* __recombination__: Mean of a Poisson distributions for number of crossing overs per sexual reproduction.
* __initial\_energy__: A parameter for parental care of infants. Values more than 0 indicate that newly born individuals can survive for a number of times without requiring food from the environment/other species. The consumption rate is always constant for all species at one unit per time step. Use with care. Having initial energy larger than zero can lead to infinite population growth because agents without food can reproduce and their offsprings also reproduce without food. For example, this happens when start age of reproduction is 1 and initial energy is larger than 0.
* __bottlenecks__: An vector of matrices each of which has the same size as space and contains the probability of death at each site. The vector is as long as model steps (generations) plus 1 (for time zero). Use this if you want to impose a killing event in the model at specific times and places. Otherwise, just use 0.0 for all the values in matrices.
* __mating\_scheme__: One of the following: -1, 0, 1. Only used in sexual reproduction. Zero means randomly paired couple are likely to have as many offsprings as any other pair. Their phenotypes does not affect their reproduction success. -1 is disassortative mating, where the more different the phenotypes of the pair, the more children they have. 1 is the opposite, assortative mating.
* __abiotic\_variance__: Variance of a normal distribution used in determining the phenotypic distance of agents to the optimal environmental phenotypes. The larger the variance, the less important is the distance.
* __biotic\_variance__: Same as _abiotic\_variance_ but for determining the biotic phenotypic distance between two individuals (used in any kind of interaction). Larger values mean that all pairs are equally likely to interact, regardless of their phenotype difference.

### 2. Model parameters

This dictionary should be named "model_parameters".

* __species__: A vector containing the name of the species _dictionaries_.
* __generations__: An _integer_ for the number of steps the model will run.
* __space__: A _tuple of two integers_ (e.g. `(3,2)`) that determines the size space. If you do not want any spatial structure, use `(1,1)`.
* __metric__: Either "chebyshev" or "euclidian". Determines how many neighbors a space site has. "chebyshev" metric means that the r-neighborhood of a position are all positions within the hypercube having side length of 2*floor(r) and being centered in the origin position. "euclidean" metric means that the r-neighborhood of a position are all positions whose cartesian indices have Euclidean distance ≤ r from the cartesian index of the given position.
* __periodic__: _Boolean__ (true or false) to determine whether boundaries of the space are connected or not.
* __resources__: A _vector of matrices_ where each matrix is the size of space and determines the amount of resources per site. The vector is as long as model steps (generations) plus 1 (for time zero).
* __interactions__: A species-species interaction _matrix of floating numbers_ determining how individuals from different species interact. Each value  is strength of interaction (between 0 and 1). Sign (+/-) is the direction of interaction where positive means similar individuals interact more strongly and negative is dissimilar ones tend to interact more.
* __food\_sources__: A species-species food _matrix of floating numbers_ determining what each species feeds on (consumption rate). Non-zero diagonal means the food resource is from the environment. Off-diagonals mean an species (in the rows) feeds on another species (in the columns). Numbers can be zero or any positive number. The magnitude of the number determines how much energy the predator gains by eating the prey. All species use 1 unit of energy per time step.
* __seed__: Either an _integer_ or `nothing` for a random number.

## Simulation outline

The simulations are fully agent-based, meaning that agents do not receive any model-level knowledge for what happens to them. The following steps happen in order to agents that are activated in a random order.

1. Grow one step older.
2. Migrate.
3. Burn energy.
4. Eat from the environment if the species can.
5. Interact with other individuals.
6. Reproduce.
7. Survive. Agents die if they are too old, or do not have enough energy, or by chance given its fitness.
8. Kill agents will a probability given in `bottlenecks`.

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

$ | 0.5 - \text{cdf}(\text{Normal}(\text{ph1},  1), \text{ph2})| $

, where cdf is cumulative density function, Normal is a normal distribution with mean ph1 (phenotypic value 1) and variance 1, and ph2 is phenotypic value 2.

If the two individuals are not predator-prey, and they are from the same species but different sexes and they are both in reproduction age, then they are marked to reproduce in the next stage. Otherwise, if they are from different species and the two species interact with each other (the corresponding values in the _interactions_ matrix is non-zero) they interact with each other. If the individual 1 to 2 has a positive value in the _interactions_ matrix, then it increases the fitness of individual 2. If the element is zero, it does not increase or decrease the fitness of individual 2. And it the element is negative, it decreases the fitness of individual 2. Similarly, individual 2 can affect the fitness of individual 1. The two interactions do not need to be symmetric. The magnitude of change in the fitness due to interaction is determined by the phenotypic distance between the two individuals and interaction value in the _interactions_ matrix.

### Reproduction

If an individual is haploid and is in reproduction age, it reproduces an identical offspring to itself except that the expression array, pleiotropy matrix, and the epistasis matrix of the offspring will mutate given the probabilities and magnitudes in the _mutation probabilities_ and _mutation magnitudes_ matrices. The number of offsprings is a random number from a Poisson distribution with a mean equal to the species' growth rate times the fitness of the individual.

If two diploid individuals are mated to reproduce, their reproduction success is proportional to their phenotypic distance (depending on `mating\_scheme`). Number of offsprings is a random number from a Poisson distribution with mean equal to the reproduction success of the two individuals times the growth rate of the species.

Each offspring of diploid individuals inherits a gamete from each parent. A gamete is a half of the expression array, pleiotropy matrix, and epistasis matrix. If recombination is allowed (> 0), then the gametes undergo crossing over. The number of crossing overs is a random number from a Poisson distribution with mean equal to the _recombination_ parameter.

#### Survival

At each step, the agent may die due to several factors. First, it dies if it has negative energy (has not been able to eat enough). Second, it dies if it is too old (> max age). Third, it dies with a probability negatively correlated to its fitness. This probability is adjusted to the _selection coefficient_. If the selection coefficient is zero, then the individual does not die. If it is one, the fitness determines survival 100%.

#### Fitness

Fitness of individuals is their abiotic fitness, i.e. their abiotic phenotypic distance to the optimal phenotypes in the environment. As such, their fitness can be different at different sites.

#### Mutation

Mutation can happen at three levels: changing the expression of each gene, changing the pleiotropy matrix, and changing the epistasis interactions between genes. The probability that a mutation occurs at each of these levels is controlled by parameter _mutation probabilities_. And size of mutations when they occur are controlled by the _mutation magnitudes_ parameters, which are the variances of a normal distribution. A mutation in gene expression and epistasis matrix is carried by adding a random number from a normal distribution to the existing values. The pleiotropy matrix is mutated by switching 0s and 1s with the probability given in the third element of _mutation magnitudes_.
