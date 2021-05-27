# 1. Make migration a personal trait affected by genes.

* [x] Users should define migration trait and specify its genes with pleiotropy.
* [x] We cannot have an optimal phenotype for migration to be selected, right?
* [x] covMat, however, should include migration trait.
* [x] users should provide the numbers and the type of distribution.
* [x] remove migration rates? replace it with Physical distances
* [ ] Migration cost? a function of fitness or the migration probability itself.
* [ ] TODO: what to do when migration trait is > 1?

# 2. Variable fitness peaks per site and per species.

## Plan 1

* TODO: Our epistasis right now has a problem: to have the following to the individual phenotypes:
  * Fmat = agent.pleiotropyMat * (agent.epistasisMat * agent.q)
  * But here we are doing unnecessary computations by multiplying complete `epistasisMat` to q and then multiplying it by `pleiotropyMat`. This can be summarized by only multiplying `pleiotropyMat` and `q` but this time th `pleiotropyMat` is not only zeros and ones, but zeros and values of `epistasisMat`. 
* calculate optPhenotypes per individual.
* growthRates should be connected to genotype? It should be equal to fitness. Population growth rate is the sum of all individual growth rates.
Keep the fitness function and population growth function as they are, separate. At each step, update the fitness peaks for each species, given the other species.

# 3. Make recombination optional and give it a prob
