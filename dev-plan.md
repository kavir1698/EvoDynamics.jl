# 1. Make migration a personal trait affected by genes.

* Users should define migration trait and specify its genes with pleiotropy.
* We cannot have an optimal phenotype for migration to be selected, right?
* covMat, however, should include migration trait
* E: users should provide the numbers and the type of distribution.
* growthRates should be connected to genotype?
* remove migration rates?

# 2. Variable fitness peaks per site and per species.

## Plan 1

Keep the fitness function and population growth function as they are, separate. At each step, update the fitness peaks for each species, given the other species.

