var documenterSearchIndex = {"docs":
[{"location":"example2/","page":"Predator prey","title":"Predator prey","text":"EditURL = \"https://github.com/kavir1698/EvoDynamics.jl/blob/master/examples/example2.jl\"","category":"page"},{"location":"example2/#Predator-prey","page":"Predator prey","title":"Predator prey","text":"","category":"section"},{"location":"example2/","page":"Predator prey","title":"Predator prey","text":"Here we create a predator prey model of two species, one haploid and one diploid.","category":"page"},{"location":"example2/","page":"Predator prey","title":"Predator prey","text":"using EvoDynamics","category":"page"},{"location":"example2/","page":"Predator prey","title":"Predator prey","text":"Model parameters are in a YAML file as follows:","category":"page"},{"location":"example2/","page":"Predator prey","title":"Predator prey","text":"species:\n  1:\n    name: a\n    number of genes: 7\n    number of phenotypes: 4\n    abiotic phenotypes: [1]\n    biotic phenotypes: [2, 3]\n    migration phenotype: 4  # can be 0 for no migration\n    migration threshold: 3.5  # phenotypic threshold after which migration is possible\n    vision radius: 1  # the radius of neighboring sites that the agent can see\n    check fraction: 0.5  # the fraction of the observable sites that the agent checks\n    ploidy: 2\n    # A random epistasis matrix where the diagonal is 1.0, meaning that each locus affects itself 100%.\n    # ngenes x (ngenes x ploidy)\n    epistasis matrix: [1.0, 0.34, 0.05, -0.12, 0.25, -0.47, 0.37, -0.32, 0.02, 0.41, 0.03, 0.0, 0.46, -0.22, 0.43, 1.0, 0.27, 0.33, 0.25, 0.19, 0.09, 0.15, -0.16, -0.18, 0.25, 0.44, -0.1, 0.21, -0.41, -0.19, 1.0, -0.48, -0.14, -0.24, -0.11, 0.4, 0.33, -0.16, 0.14, -0.34, 0.14, 0.46, 0.38, -0.36, 0.04, 1.0, 0.49, 0.41, 0.4, -0.24, 0.48, -0.4, -0.36, -0.42, -0.22, -0.01, 0.48, 0.38, 0.01, -0.4, 1.0, 0.08, 0.42, -0.21, -0.42, 0.01, 0.28, 0.37, -0.26, -0.35, -0.43, -0.28, -0.14, -0.48, 0.28, 1.0, 0.45, 0.5, 0.39, 0.04, -0.18, -0.04, 0.13, -0.11, -0.1, 0.24, 0.3, -0.22, -0.34, 0.11, 1.0, 0.22, 0.2, 0.07, 0.09, 0.43, -0.5, 0.25, -0.08, -0.22, -0.28, -0.36, -0.49, 0.03, -0.01, 1.0, -0.11, 0.2, -0.2, -0.25, -0.41, -0.03, -0.09, 0.12, 0.43, -0.24, 0.45, 0.15, -0.47, -0.33, 1.0, -0.37, 0.46, 0.21, -0.31, 0.18, -0.5, 0.12, -0.13, -0.07, -0.14, 0.49, 0.07, 0.48, 0.46, 1.0, -0.48, 0.19, 0.0, -0.38, 0.41, -0.12, 0.2, -0.12, 0.26, 0.04, 0.5, -0.49, -0.06, -0.33, 1.0, 0.29, -0.15, -0.4, 0.44, -0.39, -0.02, -0.49, -0.13, 0.41, 0.44, 0.07, 0.22, 0.49, -0.21, 1.0, 0.29, -0.28, -0.21, 0.21, 0.25, -0.37, -0.44, -0.19, -0.18, 0.5, -0.3, 0.05, 0.41, -0.02, 1.0, 0.05, -0.12, 0.26, -0.39, 0.27, -0.17, 0.13, -0.2, -0.07, 0.31, -0.42, -0.46, 0.06, 0.17, 1.0]\n    # nphenotype x (ngenes x ploidy)\n    pleiotropy matrix: [1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0]\n    growth rate: 0.8  # max number of offsprings per mating mean of a Poisson\n    expression array: [0.28878032859775615, 0.4629421231828499, 0.26092147517051467, 0.952859489607121, 0.9638502824424, 0.05038142018016245, 0.05930756376654234, 0.033459292878885716, 0.32421526342800044, 0.9029235877297073, 0.7670060809312949, 0.12766808941531993, 0.8656895869985795, 0.342191940658253]  # ngenes x ploidy long\n    selection coefficient: 0.5\n    mutation probabilities: [0.9, 0.0, 0.0]  # for gene expression array, pleiotropy matrix and epistasis matrix, respectively\n    mutation magnitudes: [0.05, 0.0, 0.01]  # same as above\n    N: [1000, 0, 0, 0]  # number of individuals per site at time 0\n    environmental noise: 0.01  # variance of a normal distribution with mean 0\n    # each row is the optimal phenotypes for each site for all abiotic traits. There are as many element as number of sites times number of abiotic traits. The first N elements are for the first abiotic trait, where N is the number of sites, and so on.\n    optimal phenotype values:\n      - [2.8184972630154848, 1.2669190502061671, 1.7947315210311097, 0.5627451656026976]\n      - [2.3186137078797455, 0.7146284070374198, 1.7331797945458374, 1.6353360768904688]\n    # Optimal phenotype indices for each generation, including generation zero.\n    optimal phenotypes: [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2]\n    age: 5  # max age\n    recombination: 1 # Mean of a Poisson distributions for number of crossing overs\n    initial energy: 0 # A parameter for parental care of infants. Values more than 0 indicate that newly born individuals can survive for a number of times without requiring food from the environment/other species. The consumption rate (i.e. how many generations this initial energy suffices) is determined by the sum of the corresponding rows in \"food sources\"\n\n  2:\n    name: b\n    number of genes: 8\n    number of phenotypes: 5\n    abiotic phenotypes: [1,2]\n    biotic phenotypes: [3, 4]\n    migration phenotype: 5\n    migration threshold: 3.4\n    vision radius: 1\n    check fraction: 0.5\n    ploidy: 1\n    epistasis matrix: [1.0, 0.28, 0.0, 0.0, 0.31, -0.19, -0.43, 0.26, -0.01, 1.0, 0.28, -0.42, 0.47, 0.29, 0.38, 0.27, -0.01, -0.21, 1.0, 0.05, -0.42, -0.33, -0.06, -0.44, -0.34, -0.11, 0.1, 1.0, -0.34, -0.49, 0.39, -0.08, -0.42, 0.12, 0.09, -0.11, 1.0, 0.17, 0.21, 0.47, 0.18, -0.46, 0.3, -0.07, -0.4, 1.0, -0.5, -0.27, -0.09, 0.21, 0.1, -0.27, 0.16, -0.1, 1.0, -0.27, 0.13, -0.39, 0.17, 0.43, 0.11, -0.28, -0.08, 1.0]\n    pleiotropy matrix: [1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1]\n    covariance matrix: [0.11035, 0.2686, 0.6074, 0.6914, 0.8604, 0.2686, 0.6807, 0.1738, 0.619, 0.1465, 0.6074, 0.1738, 0.6494, 0.7646, 0.75, 0.6914, 0.619, 0.7646, 0.4658, 0.549, 0.8604, 0.1465, 0.75, 0.549, 0.6465]\n    growth rate: 1.2\n    expression array: [0.24923147816626035, 0.7155732641738595, 0.9655184311211502, 0.8638149724268844, 0.5075272565823061, 0.9189652626508431, 0.7897640036022151, 0.17091233765481717]\n    selection coefficient: 0.5\n    mutation probabilities: [0.9, 0.0, 0.0]\n    mutation magnitudes: [0.05, 0.0, 0.01]\n    N: [1000, 0, 0, 0]\n    environmental noise: 0.01\n    optimal phenotype values:\n      - [0.7614758101208934, 2.2091361414343313, 0.7920974352892358, 2.587205421882114, 2.3911353663016586, 1.7858661540288683, 0.7630236717263839, 2.311211631439866]\n      - [0.6437641305315445, 0.9954545836033715, 2.469792530385348, 1.6158867433451882, 2.097629262577973, 1.1314248848669466, 1.490299526620522, 0.2566056477862022]\n    optimal phenotypes: [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2]\n    age: 3\n    recombination: 1 # Mean of a Poisson distributions for number of crossing overs\n    initial energy: 0\n\nmodel:\n  generations: 14  # number of simulation steps\n  space: [2, 2]\n  metric: chebyshev  # how many neighbors a space site has. \"chebyshev\" metric means that the r-neighborhood of a position are all positions within the hypercube having side length of 2*floor(r) and being centered in the origin position. \"euclidean\" metric means that the r-neighborhood of a position are all positions whose cartesian indices have Euclidean distance ≤ r from the cartesian index of the given position.\n  periodic: false  # whether boundaries of the space are connected\n  resources: [200, 158, 183, 190]  # available resources per site\n  interactions: [0.1, 0.0, 0.0, -1.0]  # How individuals from different species interact. value  is strength of interaction (between 0 and 1). Sign is the direction of interaction where positive means similar individuals interact more strongly and negative is dissimilar ones tend to interact more.\n  food sources: [1.0, 0.5, 0.0, 0.0]  # What each species feeds on (consumption rate). Has priority over interactions. Non-zero diagonal means the food resource is from the environment. It will be read from rows (species order) to columns (species order).\n  seed: 2","category":"page"},{"location":"example2/","page":"Predator prey","title":"Predator prey","text":"param_file = \"../../examples/paramfile2.yml\"\nagentdata, modeldata, model = runmodel(param_file);\n\nmodeldata","category":"page"},{"location":"example2/","page":"Predator prey","title":"Predator prey","text":"","category":"page"},{"location":"example2/","page":"Predator prey","title":"Predator prey","text":"This page was generated using Literate.jl.","category":"page"},{"location":"model_description/#Model-description","page":"Model description","title":"Model description","text":"","category":"section"},{"location":"model_description/#Parameters","page":"Model description","title":"Parameters","text":"","category":"section"},{"location":"model_description/#Species-specific-parameters","page":"Model description","title":"Species specific parameters","text":"","category":"section"},{"location":"model_description/","page":"Model description","title":"Model description","text":"Each species should have the following parameters. The order that you write these parameters does not matter.","category":"page"},{"location":"model_description/","page":"Model description","title":"Model description","text":"name: a name for the species.\nnumber of genes: An integer for number of genes that the species has.\nploidy: Either 1 for haploid or 2 for diploid genomes. Diploids can have recombination.\nnumber of phenotypes: An integer for the number of phenotypes that the species has.\nabiotic phenotypes: An array of integers (e.g. \"[1,2]\") specifying abiotic phenotypes among all phenotypes. Abiotic phenotypes determine how the species interacts with the environment.\nbiotic phenotypes: An array of integers (e.g. \"[3]\") specifying biotic phenotypes among all phenotypes. Biotic phenotypes determine how the species interacts with other individuals from the same or different species.\nmigration phenotype: An integer specifying the phenotype that determines migration trait. If the species does not migrate, put 0.\nmigration threshold: The phenotypic value of the migration phenotype after which an individual can migrate.\nvision radius: A number determining the radius of neighboring sites that the agent can see before migration.\ncheck fraction: A number between 0 and 1 showing the fraction of the visible sites to the agent that it can check and decide whether to migrate to.\nepistasis matrix: An epistasis matrix is of size l times l, where l is the product of number of genes and ploidy. Epistasis matrix specifies the direction (positive or negative) and size of effect of one locus on other loci. For example, if at row 1 and column 2 is a value 0.2, it means that locus 1 affects locus 2 by increasing the effect of locus 2 (because its positive) with 20% of the effect of locus 1.\npleiotropy matrix: A binary matrix (0s and 1s) with size number of phenotypes times l. The pleiotropy matrix specifies the phenotypes that each locus affects.\nexpression array: A vector of size l that represent the expression amount of each locus determining its effect size.\ngrowth rate: Mean of a Poisson distribution for number of offsprings per reproduction. This number is the maximum mean when fitness of a haploid individual is 1, or the distance between the biotic phenotypes of two diploid individuals is 0.\nselection coefficient: A number between 0 and 1 that determines the importance of fitness. 0 would be a model without selection.\nmutation probabilities: A vector of three numbers each of which specifies the probability for a different type of mutations: mutation probability of the expression array, pleiotropy matrix, and epistasis matrix, respectively.\nmutation magnitudes: A vector of numbers with the same size as mutation probabilities that determines the magnitude of mutation for each of the three categories. Specifically, the numbers are the variances of normal distributions with mean 0 for expression array and epistasis matrices, and probability of changing a 0 and 1 in in the pleiotropy matrix.\nN: A vector of integers for the initial number of individuals at each site.\nenvironmental noise: A number for the variance of a normal distribution with mean 0 that will be added to the phenotypes.\noptimal phenotype values: One or more vectors of numbers. Each vector should start with a dash and one indentation level. Each vector is the optimal phenotypes for each site for all abiotic traits. There are as many element as number of sites times number of abiotic traits. The first N elements are for the first abiotic trait, where N is the number of sites, and so on.\noptimal phenotypes: A vector of integers for optimal phenotype indices for each generation, including generation zero.\nage: An integer for maximum age of individuals of this species.\nreproduction start age: The age at which individuals can reproduce.\nreproduction end age: The age after which individuals cannot reproduce.\nrecombination: Mean of a Poisson distributions for number of crossing overs per sexual reproduction.\ninitial energy: A parameter for parental care of infants. Values more than 0 indicate that newly born individuals can survive for a number of times without requiring food from the environment/other species. The consumption rate (i.e. how many generations this initial energy suffices) is determined by the sum of the corresponding rows in \"food sources\" model parameter.\nbottleneck function: Forcefully kill specific agents. Path to a file that has a function named \"bottleneck\". The function accepts two arguments: agent and model. It returns true or false for death or survival of the individual, respectively.\nbottleneck times: A boolean matrix (zeros and ones) turned into a vector where rows are sites and columns are time steps.  Determines time steps at which the bottleneck function should be activated.","category":"page"},{"location":"model_description/#Model-parameters","page":"Model description","title":"Model parameters","text":"","category":"section"},{"location":"model_description/","page":"Model description","title":"Model description","text":"generations: An integer for the number of steps the model will run.\nspace: A vector of two integers that determine the size of a grid for space.\nmetric: Either \"chebyshev\" or \"euclidian\". Determines how many neighbors a space site has. \"chebyshev\" metric means that the r-neighborhood of a position are all positions within the hypercube having side length of 2*floor(r) and being centered in the origin position. \"euclidean\" metric means that the r-neighborhood of a position are all positions whose cartesian indices have Euclidean distance ≤ r from the cartesian index of the given position.\nperiodic: _Boolean__ (true or false) to determine whether boundaries of the space are connected or not.\nresources: A vector of integers determining available resources (e.g. vegetation) per site per time step.\ninteractions: A species-species interaction matrix of numbers determining how individuals from different species interact. Each value  is strength of interaction (between 0 and 1). Sign (+/-) is the direction of interaction where positive means similar individuals interact more strongly and negative is dissimilar ones tend to interact more.\nfood sources: A species-species food matrix of numbers determining what each species feeds on (consumption rate). Non-zero diagonal means the food resource is from the environment. Non-diagonals mean an species (in the rows) feeds on another species (in the columns). Numbers can be zero or any positive number. The magnitude of the number determines how many generations can an individual live off of given one unit of the food source. For example, if a diagonal is 2, it means that the species will eat one unit of the environmental resources and that is enough for it to live two steps.\nseed: Either an integer or Null for random number generator seed.","category":"page"},{"location":"model_description/#Simulation-outline","page":"Model description","title":"Simulation outline","text":"","category":"section"},{"location":"model_description/","page":"Model description","title":"Model description","text":"The simulations are fully agent-based, meaning that agents do not receive any model-level knowledge for what happens to them. The following steps happen in order to agents that are activated in a random order.","category":"page"},{"location":"model_description/","page":"Model description","title":"Model description","text":"Grow one step older.\nMigrate.\nBurn energy.\nEat from the environment if the species can.\nInteract with other individuals.\nIf meeting another individual of the same species but with different sex, try reproduction.\nIf haploid, reproduce.\nSurvive. Agents die if they are too old, or do not have enough energy, or by chance given its fitness.\nGo through the bottleneck function if it should be activated at the current time.","category":"page"},{"location":"model_description/#Migration","page":"Model description","title":"Migration","text":"","category":"section"},{"location":"model_description/","page":"Model description","title":"Model description","text":"If a species have the capability to migrate (i.e. migration phenotype is not 0), then at each time step, it is checked whether the phenotypic value of the migration phenotype is above migration threshold. If it is, then the agent checks a random check fraction of the sites neighboring it and moves to the most suitable site. Suitability is calculated by comparing the agent's abiotic phenotype and the optimal abiotic phenotypic value for each site. If the checked neighbors have worse conditions than the current site, the agent does not move.","category":"page"},{"location":"model_description/#Burn-energy","page":"Model description","title":"Burn energy","text":"","category":"section"},{"location":"model_description/","page":"Model description","title":"Model description","text":"Agent's energy is reduced by one unit.","category":"page"},{"location":"model_description/#Eat","page":"Model description","title":"Eat","text":"","category":"section"},{"location":"model_description/","page":"Model description","title":"Model description","text":"If the agent is able to eat from the environment (the diagonal of food sources at the corresponding row and column is non-zero) and there is any environmental resources left at the agent's site, then its energy level is boosted by as much as the value in the corresponding element at food sources.","category":"page"},{"location":"model_description/#Interacting-with-other-individuals","page":"Model description","title":"Interacting with other individuals","text":"","category":"section"},{"location":"model_description/","page":"Model description","title":"Model description","text":"At each time, the agent interacts with as many agents as there are species in the model. If there are fewer agents in the site than number of species, the agent interacts with all of them. Otherwise, it interacts with at most one randomly picked individual from each species. If most of the individuals at the site are from the first species, then it is more likely that the agent interacts with an individual from the first species and with no individual from the second species. This set up allows interactions between species be dependent on the population size of the species at each site.","category":"page"},{"location":"model_description/","page":"Model description","title":"Model description","text":"For each pair of interacting individuals, we first check whether one individual feeds on the other one. If so, the hunt is successful with a probability proportional to the average phenotypic distance between the biotic phenotypes of the two individuals.","category":"page"},{"location":"model_description/","page":"Model description","title":"Model description","text":"Phenotypic distance between two species depends on the sign of the corresponding element in interactions matrix. If it is negative, then the hunt is more successful if the average phenotype of the two individuals is more different.","category":"page"},{"location":"model_description/","page":"Model description","title":"Model description","text":"The phenotypic distance between two individuals is the average distance between all pairs of biotic phenotypes between the two individuals. To calculate the distance between two phenotypic values, we use the following formula:  ","category":"page"},{"location":"model_description/","page":"Model description","title":"Model description","text":"$ | 0.5 - \\text{cdf}(\\text{Normal}(\\text{ph1},  1), \\text{ph2})| $","category":"page"},{"location":"model_description/","page":"Model description","title":"Model description","text":", where cdf is cumulative density function, Normal is a normal distribution with mean ph1 (phenotypic value 1) and variance 1, and ph2 is phenotypic value 2.","category":"page"},{"location":"model_description/","page":"Model description","title":"Model description","text":"If the two individuals are not predator-prey, and they are from the same species but different sexes and they are both in reproduction age, then they try to reproduce. Otherwise, if they are from different species and the two species interact with each other (the corresponding values in the interactions matrix is non-zero) they interact with each other. If the individual 1 to 2 has a positive value in the interactions matrix, then it increases the fitness of individual 2.  If the element is zero, it does not increase or decrease the fitness of individual 2. And it the element is negative, it decreases the fitness of individual 2. Similarly, individual 2 can affect the fitness of individual 1. The two interactions do not need to be symmetric. The magnitude of change in the fitness due to interaction is determined by the phenotypic distance between the two individuals and interaction value in the interactions matrix.","category":"page"},{"location":"model_description/#Reproduction","page":"Model description","title":"Reproduction","text":"","category":"section"},{"location":"model_description/","page":"Model description","title":"Model description","text":"If an individual is haploid and is in reproduction age, it reproduces an identical offspring to itself except that the expression array, pleiotropy matrix, and the epistasis matrix of the offspring will mutate given the probabilities and magnitudes in the mutation probabilities and mutation magnitudes matrices. The number of offsprings is a random number from a Poisson distribution with a mean equal to the species' growth rate times the fitness of the individual.","category":"page"},{"location":"model_description/","page":"Model description","title":"Model description","text":"If two diploid individuals are mated to reproduce, their reproduction success is proportional to their phenotypic similarity. Number of offsprings is a random number from a Poisson distribution with mean equal to the reproduction success of the two individuals times the growth rate of the species.","category":"page"},{"location":"model_description/","page":"Model description","title":"Model description","text":"Each offspring of diploid individuals inherits a gamete from each parent. A gamete is a half of the expression array, pleiotropy matrix, and epistasis matrix. If recombination is allowed (> 0), then the gametes undergo crossing over. The number of crossing overs is a random number from a Poisson distribution with mean equal to the recombination parameter.","category":"page"},{"location":"model_description/#Survival","page":"Model description","title":"Survival","text":"","category":"section"},{"location":"model_description/","page":"Model description","title":"Model description","text":"At each step, the agent may die due to several factors. First, it dies if it has negative energy (has not been able to eat enough). Second, it dies if it is too old (> max age). Third, it dies with a probability negatively correlated to its fitness. This probability is adjusted to the selection coefficient. If the selection coefficient is zero, then the individual does not die. If it is one, the fitness determines survival 100%.","category":"page"},{"location":"model_description/#Fitness","page":"Model description","title":"Fitness","text":"","category":"section"},{"location":"model_description/","page":"Model description","title":"Model description","text":"Fitness of individuals is their abiotic fitness, i.e. their abiotic phenotypic distance to the optimal phenotypes in the environment. As such, their fitness can be different at different sites.","category":"page"},{"location":"model_description/#Mutation","page":"Model description","title":"Mutation","text":"","category":"section"},{"location":"model_description/","page":"Model description","title":"Model description","text":"Mutation can happen at three levels: changing the expression of each gene, changing the pleiotropy matrix, and changing the epistasis interactions between genes. The probability that a mutation occurs at each of these levels is controlled by parameter mutation probabilities. And size of mutations when they occur are controlled by the mutation magnitudes parameters, which are the variances of a normal distribution. A mutation in gene expression and epistasis matrix is carried by adding a random number from a normal distribution to the existing values. The pleiotropy matrix is mutated by switching 0s and 1s with the probability given in the third element of mutation magnitudes.","category":"page"},{"location":"example1/","page":"Simple Wright-Fisher","title":"Simple Wright-Fisher","text":"EditURL = \"https://github.com/kavir1698/EvoDynamics.jl/blob/master/examples/example1.jl\"","category":"page"},{"location":"example1/#Simple-Wright-Fisher","page":"Simple Wright-Fisher","title":"Simple Wright-Fisher","text":"","category":"section"},{"location":"example1/","page":"Simple Wright-Fisher","title":"Simple Wright-Fisher","text":"We can create and run simple Wright-Fisher simulations with EvoDynamics.jl. To that end, we define a single haploid species, in an unstructured space, with two single genes affecting biotic and abiotic traits, respectively.","category":"page"},{"location":"example1/","page":"Simple Wright-Fisher","title":"Simple Wright-Fisher","text":"using EvoDynamics","category":"page"},{"location":"example1/","page":"Simple Wright-Fisher","title":"Simple Wright-Fisher","text":"A simple one-species model with no spatial structure. Model parameters are in a YAML file as follows:","category":"page"},{"location":"example1/","page":"Simple Wright-Fisher","title":"Simple Wright-Fisher","text":"species:\n  1:\n    name: a\n    number of genes: 2\n    number of phenotypes: 2\n    abiotic phenotypes: [1]\n    biotic phenotypes: [2]\n    migration phenotype: 0\n    migration threshold: 3.5\n    vision radius: 0\n    check fraction: 0\n    ploidy: 1\n    epistasis matrix: [1.0, 0.0, 0.0, 1.0]\n    pleiotropy matrix: [1, 0, 0, 1]\n    growth rate: 1.0\n    expression array: [0.28, 0.46]\n    selection coefficient: 0.5\n    mutation probabilities: [0.9, 0.0, 0.0]\n    mutation magnitudes: [0.05, 0.0, 0.0]\n    N: [100]\n    environmental noise: 0.01\n    optimal phenotype values:\n      - [1.76]\n    optimal phenotypes: [1, 1, 1, 1, 1, 1]\n    age: 2\n    recombination: 0\n    initial energy: 0\n\nmodel:\n  generations: 5\n  space: [1,1]\n  metric: chebyshev\n  periodic: false\n  resources: [200]\n  interactions: [-0.1]\n  food sources: [1.0]\n  seed: Null","category":"page"},{"location":"example1/","page":"Simple Wright-Fisher","title":"Simple Wright-Fisher","text":"param_file = \"../../examples/paramfile1.yml\" #hide\nagentdata, modeldata, model = runmodel(param_file);\n\nmodeldata","category":"page"},{"location":"example1/","page":"Simple Wright-Fisher","title":"Simple Wright-Fisher","text":"","category":"page"},{"location":"example1/","page":"Simple Wright-Fisher","title":"Simple Wright-Fisher","text":"This page was generated using Literate.jl.","category":"page"},{"location":"#EvoDynamics.jl-Documentation","page":"Introduction","title":"EvoDynamics.jl Documentation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"EvoDynamics.jl tries to bridge the gap in studying biological systems at small and large scales. Some studies only focus on single genes affecting single phenotypes, some studies only analyze gene interactions, some focus on populations, and some on species interactions. EvoDynamics.jl is a framework to study the effect of interactions between all these levels. It includes explicit pleiotropy, epistasis, selection acting on multiple phenotypes, different phenotypes affecting fitness at different amounts, arbitrary spatial structure, migration, and interacting species.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Figure below shows different biological levels controlled by the model.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"(Image: Fig. 1. __Model structure.__)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"See Tutorial for running the model, and Model description for a description of model parameters and simulation outline.","category":"page"},{"location":"#Features","page":"Introduction","title":"Features","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Possibility to model complex food webs with various asymmetrical interactions.\nIndividuals interact given their phenotypes.\nConnecting genome structure to phenotypes to populations.\nSpace is a grid on which individuals can migrate.\nPossibility to model a complex environment with various amounts of resources at each site.\nAllows time varying optimal phenotypic value per site.\nPossibility of killing specific individuals at certain times and sites.\nCan model both haploid and diploid species.\nIncludes mutation and recombination.\nEasy data collection.\nRuns replicates and collects data in a table automatically.\nCan run replicated in parallel.","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/#Installation","page":"Tutorial","title":"Installation","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Install using the following command in a Julia REPL.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"]add EvoDynamics","category":"page"},{"location":"tutorial/#Basic-usage","page":"Tutorial","title":"Basic usage","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Parameters of a model should be put in a YAML file with the structure below. Note that spaces and indentations are meaningful in YAML. Indentations should be spaces not tabs.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"See Simple Wright-Fisher and Predator prey for complete examples of parameter files.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"species:\n  1:\n    name: a\n    parameter 1: ...\n    parameter 2: ...\n    ...\n  2:\n    name: b\n    parameter 1: ... \n    parameter 2: ... \n    ...\nmodel:\n  model parameter 1: ...\n  model parameter 2: ...\n  ...","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The file has two main levels: species and model. species stores species specific parameters as many different species as you want.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Since we cannot write a matrix in a YAML file, any parameter that is a matrix should be converted to a vector. In Julia, you can do this by vec(yourmatrix).","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"First, define your model parameters in a YAML file (here, we call it parameters.yml). Simple Wright-Fisher and Predator prey have examples of initiation parameters. See Model description for a description of each parameter.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"We can the use the runmodel function to create a model from these parameters and run the simulation.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"runmodel","category":"page"},{"location":"tutorial/#EvoDynamics.runmodel","page":"Tutorial","title":"EvoDynamics.runmodel","text":"runmodel(param_file::AbstractString; kwargs)\n\nCreates and runs a model given parameters. Returns a DataFrame of collected data, which are specified by kwargs.\n\nKeywords\n\nadata=[] agent data to be collected. Either agent fields or functions that accept an agent as input can be put in the array. To aggregate collected data, provide tuples inside the array. For example, to collect mean and median fitness of individuals which is in field W, your array will be [(:W,mean), (:W,median)].\nmdata=[meanfitnessper_species] model data to be collected. By default, collects mean population fitness per species. Each row of the output DataFrame corresponds to all agents and each column is the value function applied to a field. The functions in the array are applied to the model object.\nwhen=nothing The generations from which data are collected. By default collect at all steps.\nreplicates::Int = 0 Number of replicates per simulation.\nparallel::Bool = false Whether to run replicates in parallel. If true, you should add processors to your julia session (e.g. by addprocs(n)) and define your parameters and EvoDynamics on all workers. To do that, add @everywhere before them. For example, @everywhere EvoDynamics.\nseeds = optionally, provide an array of integers as seeds for each replicate.\n\n\n\n\n\n","category":"function"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using EvoDynamics\nagentdata, modeldata, model = runmodel(\"parameters.yml\")","category":"page"},{"location":"tutorial/#Creating-simulation-parameter-files","page":"Tutorial","title":"Creating simulation parameter files","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"EvoDynamics.jl reads simulation parameters (Model description) from a human-readable YAML file. This file can be populated manually using a text editor or from within a Julia session.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To create parameters from within a Julia session and write them to a YAML file, you can follow the example below.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To define species parameters, create a dictionary whose keys are numbers starting from 1 for the number of species you have, and values are dictionaries themselves with the species parameter names and values. In the example below, for brevity, I only add a few parameters for each species (see Model description for the complete list of parameters).","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"species_params = Dict(\n  1=>Dict(\"name\" => \"a\", \"number of genes\" => 2, \"number of phenotypes\" => 2, \"abiotic phenotypes\" => [1]),\n  2=>Dict(\"name\" => \"b\", \"number of genes\" => 2, \"number of phenotypes\" => 2, \"abiotic phenotypes\" => [1])\n)\n","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To define model parameters, create a new dictionary whose keys are the parameter names and values the parameter values. Here, again I only define a few parameters for brevity.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"\nmodel_params = Dict(\"generations\"=> 100, \"space\"=> [6,10], \"food sources\" => [1.0, 0.7, 0.0, 0.0])\n","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Finally, create a dictionary mixing the two dictionaries before with keys \"species\" and \"model\". ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"data = Dict(\"species\"=> species_params, \"model\"=> model_params)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"This dictionary can be written as a YAML file in the correct format.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using YAML\nf = \"params.yml\"\nYAML.write_file(f, data)","category":"page"},{"location":"tutorial/#Collecting-data","page":"Tutorial","title":"Collecting data","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The interface to the model is from the runmodel function (see Tutorial).","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"EvoDynamics.jl uses Agents.jl underneath. See Agents.jl's documentation for details about writing functions to collect any data during simulations. Here, we explain the specific implementation of the model.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"There are two main objects from which you can collect data: and agent object of type AbstractAgent and a model object of type ABM. Both of these types are defined the Agents.jl package.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Agent object has the following fields that are individual specific: id, pos, species, epistasisMat (epistasis matrix), pleiotropyMat (pleiotropy matrix), q (gene expression array), biotic_phenotype, abiotic_phenotype, age, sex, energy, interaction_history, and W (fitness).","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The model object has the following fields that can be accessed with the . syntax and describe properties of species or the model: ngenes, nphenotypes, growthrates, selectionCoeffs, ploidy, optvals (optimal values), optinds (optval indices per generation), mutProbs (mutation probabilities), mutMagnitudes (mutation magnitudes), N, E (environmental noise), generations, nspecies, migration_traits, vision_radius, check_fraction, migration_thresholds, step, biotic_phenotypes (indices of biotic phenotypes per species), abiotic_phenotypes, max_ages, food_sources, interactions, resources, recombination, initial_energy.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"You can collect data from agents and/or from the model object. To collect data from agents, use the adata keyword argument in the runmodel function, and to collect data from the model, use the mdata keyword. A complete description of the values these keywords take are at data collection section of the Agents.jl package.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"For example, we use the function below to count the number of individual per species:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"\"Returns the size of each species.\"\nfunction species_N(model::ABM)\n  allagents = model.agents\n  if length(allagents) == 0\n    return fill(0, model.nspecies)\n  else\n    counts = countmap([a.species for a in values(model.agents)])\n    output = fill(0, model.nspecies)\n    for (k, v) in counts\n      output[k] = v\n    end\n    return output\n  end\nend\n\nusing EvoDynamics\n\nagentdata, modeldata, model = runmodel(\"parameters.yml\", mdata=[species_N])","category":"page"}]
}
