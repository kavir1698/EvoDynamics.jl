```@meta
EditURL = "<unknown>/examples/example2.jl"
```

# Predator prey

Here we create a predator prey model of two species, one haploid and one diploid.

```@example example2
using EvoDynamics
```

Model parameters are in a YAML file as follows:

```yml
- species:
  - id: 1
    number of genes: 7
    number of phenotypes: 4
    abiotic phenotypes: [1]
    biotic phenotypes: [2, 3]
    migration phenotype: 4  # can be 0 for no migration
    migration threshold: 3.5  # phenotypic threshold after which migration is possible
    vision radius: 1  # the radius of neighboring sites that the agent can see
    check fraction: 0.5  # the fraction of the observable sites that the agent checks
    ploidy: 2
    # A random epistasis matrix where the diagonal is 1.0, meaning that each locus affects itself 100%.
    # ngenes x (ngenes x ploidy)
    epistasis matrix: [1.0, 0.34, 0.05, -0.12, 0.25, -0.47, 0.37, -0.32, 0.02, 0.41, 0.03, 0.0, 0.46, -0.22, 0.43, 1.0, 0.27, 0.33, 0.25, 0.19, 0.09, 0.15, -0.16, -0.18, 0.25, 0.44, -0.1, 0.21, -0.41, -0.19, 1.0, -0.48, -0.14, -0.24, -0.11, 0.4, 0.33, -0.16, 0.14, -0.34, 0.14, 0.46, 0.38, -0.36, 0.04, 1.0, 0.49, 0.41, 0.4, -0.24, 0.48, -0.4, -0.36, -0.42, -0.22, -0.01, 0.48, 0.38, 0.01, -0.4, 1.0, 0.08, 0.42, -0.21, -0.42, 0.01, 0.28, 0.37, -0.26, -0.35, -0.43, -0.28, -0.14, -0.48, 0.28, 1.0, 0.45, 0.5, 0.39, 0.04, -0.18, -0.04, 0.13, -0.11, -0.1, 0.24, 0.3, -0.22, -0.34, 0.11, 1.0, 0.22, 0.2, 0.07, 0.09, 0.43, -0.5, 0.25, -0.08, -0.22, -0.28, -0.36, -0.49, 0.03, -0.01, 1.0, -0.11, 0.2, -0.2, -0.25, -0.41, -0.03, -0.09, 0.12, 0.43, -0.24, 0.45, 0.15, -0.47, -0.33, 1.0, -0.37, 0.46, 0.21, -0.31, 0.18, -0.5, 0.12, -0.13, -0.07, -0.14, 0.49, 0.07, 0.48, 0.46, 1.0, -0.48, 0.19, 0.0, -0.38, 0.41, -0.12, 0.2, -0.12, 0.26, 0.04, 0.5, -0.49, -0.06, -0.33, 1.0, 0.29, -0.15, -0.4, 0.44, -0.39, -0.02, -0.49, -0.13, 0.41, 0.44, 0.07, 0.22, 0.49, -0.21, 1.0, 0.29, -0.28, -0.21, 0.21, 0.25, -0.37, -0.44, -0.19, -0.18, 0.5, -0.3, 0.05, 0.41, -0.02, 1.0, 0.05, -0.12, 0.26, -0.39, 0.27, -0.17, 0.13, -0.2, -0.07, 0.31, -0.42, -0.46, 0.06, 0.17, 1.0]
    # nphenotype x (ngenes x ploidy)
    pleiotropy matrix: [1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0]
    growth rate: 0.8  # max number of offsprings per mating mean of a Poisson
    expression array: [0.28878032859775615, 0.4629421231828499, 0.26092147517051467, 0.952859489607121, 0.9638502824424, 0.05038142018016245, 0.05930756376654234, 0.033459292878885716, 0.32421526342800044, 0.9029235877297073, 0.7670060809312949, 0.12766808941531993, 0.8656895869985795, 0.342191940658253]  # ngenes x ploidy long
    selection coefficient: 0.5
    mutation probabilities: [0.9, 0.0, 0.0]  # for gene expression array, pleiotropy matrix and epistasis matrix, respectively
    mutation magnitudes: [0.05, 0.0, 0.01]  # same as above
    N: [1000, 0, 0, 0]  # number of individuals per site at time 0
    environmental noise: 0.01  # variance of a normal distribution with mean 0
    # each row is the optimal phenotypes for each site for all abiotic traits. There are as many element as number of sites times number of abiotic traits. The first N elements are for the first abiotic trait, where N is the number of sites, and so on.
    optimal phenotype values:
      - [2.8184972630154848, 1.2669190502061671, 1.7947315210311097, 0.5627451656026976]
      - [2.3186137078797455, 0.7146284070374198, 1.7331797945458374, 1.6353360768904688]
    # Optimal phenotype indices for each generation, including generation zero.
    optimal phenotypes: [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2]
    age: 5  # max age
    recombination: 1 # Mean of a Poisson distributions for number of crossing overs
    initial energy: 0 # A parameter for parental care of infants. Values more than 0 indicate that newly born individuals can survive for a number of times without requiring food from the environment/other species. The consumption rate (i.e. how many generations this initial energy suffices) is determined by the sum of the corresponding rows in "food sources"

  - id: 2
    number of genes: 8
    number of phenotypes: 5
    abiotic phenotypes: [1,2]
    biotic phenotypes: [3, 4]
    migration phenotype: 5
    migration threshold: 3.4
    vision radius: 1
    check fraction: 0.5
    ploidy: 1
    epistasis matrix: [1.0, 0.28, 0.0, 0.0, 0.31, -0.19, -0.43, 0.26, -0.01, 1.0, 0.28, -0.42, 0.47, 0.29, 0.38, 0.27, -0.01, -0.21, 1.0, 0.05, -0.42, -0.33, -0.06, -0.44, -0.34, -0.11, 0.1, 1.0, -0.34, -0.49, 0.39, -0.08, -0.42, 0.12, 0.09, -0.11, 1.0, 0.17, 0.21, 0.47, 0.18, -0.46, 0.3, -0.07, -0.4, 1.0, -0.5, -0.27, -0.09, 0.21, 0.1, -0.27, 0.16, -0.1, 1.0, -0.27, 0.13, -0.39, 0.17, 0.43, 0.11, -0.28, -0.08, 1.0]
    pleiotropy matrix: [1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1]
    covariance matrix: [0.11035, 0.2686, 0.6074, 0.6914, 0.8604, 0.2686, 0.6807, 0.1738, 0.619, 0.1465, 0.6074, 0.1738, 0.6494, 0.7646, 0.75, 0.6914, 0.619, 0.7646, 0.4658, 0.549, 0.8604, 0.1465, 0.75, 0.549, 0.6465]
    growth rate: 1.2
    expression array: [0.24923147816626035, 0.7155732641738595, 0.9655184311211502, 0.8638149724268844, 0.5075272565823061, 0.9189652626508431, 0.7897640036022151, 0.17091233765481717]
    selection coefficient: 0.5
    mutation probabilities: [0.9, 0.0, 0.0]
    mutation magnitudes: [0.05, 0.0, 0.01]
    N: [1000, 0, 0, 0]
    environmental noise: 0.01
    optimal phenotype values:
      - [0.7614758101208934, 2.2091361414343313, 0.7920974352892358, 2.587205421882114, 2.3911353663016586, 1.7858661540288683, 0.7630236717263839, 2.311211631439866]
      - [0.6437641305315445, 0.9954545836033715, 2.469792530385348, 1.6158867433451882, 2.097629262577973, 1.1314248848669466, 1.490299526620522, 0.2566056477862022]
    optimal phenotypes: [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2]
    age: 3
    recombination: 1 # Mean of a Poisson distributions for number of crossing overs
    initial energy: 0

- model:
  generations: 14  # number of simulation steps
  space: [2, 2]
  metric: chebyshev  # how many neighbors a space site has. "chebyshev" metric means that the r-neighborhood of a position are all positions within the hypercube having side length of 2*floor(r) and being centered in the origin position. "euclidean" metric means that the r-neighborhood of a position are all positions whose cartesian indices have Euclidean distance ≤ r from the cartesian index of the given position.
  periodic: false  # whether boundaries of the space are connected
  resources: [200, 158, 183, 190]  # available resources per site per time
  interactions: [0.1, 0.0, 0.0, -1.0]  # How individuals from different species interact. value  is strength of interaction (between 0 and 1). Sign is the direction of interaction where positive means similar individuals interact more strongly and negative is dissimilar ones tend to interact more.
  food sources: [1.0, 0.5, 0.0, 0.0]  # What each species feeds on (consumption rate). Has priority over interactions. Non-zero diagonal means the food resource is from the environment. It will be read from rows (species order) to columns (species order).
  seed: 2
```

```@example example2
param_file = "../../examples/paramfile2.yml"
agentdata, modeldata, model = runmodel(param_file);

modeldata
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

