# NB OBSOLETE
using Random
using LinearAlgebra: Symmetric


nphenotypes = (4, 5)
abiotic_phenotypes = ([1], [1,2])
ngenes = (7, 8)
ploidy = (2, 1)
# choose random values for epistasis matrix epistasisMat, but make sure the diagonal is 
# 1.0, meaning that each locus affects itself 100%.
epistasisMat = Tuple([Random.rand(-0.5:0.01:0.5, i, i) for i in (ngenes .* ploidy)])
for index in 1:length(epistasisMat)
  for diag in 1:size(epistasisMat[index], 1)
    epistasisMat[index][diag, diag] = 1
  end
end

generations = 5
space = (2,2)

# opt_phenotypes1.csv and opt_phenotypes2.csv are optimial phenotypes for species 1 and 2. The format is each row is the optimal phenotypes for each site for all abiotic traits. For example, species 2 has two abiotic trait and the space has four sites. opt_phenotypes2.csv has eight columns the first four of which are optimal phenotypes of the first trait and the last four of which are for the second trait.
optPhenotype_files = ("opt_phenotypes1.csv", "opt_phenotypes2.csv")
# The dictionary below assigns a row from the files above to a generation including generation 0. The keys are species numbers.
optPhenotypes = Dict(1 => [1,1,1,2,3,3], 2 => [1,1,1,2,2,2])

pleiotropyMat = (rand([true, false], nphenotypes[1], ngenes[1] * ploidy[1]), rand([true, false], nphenotypes[2], ngenes[2] * ploidy[2]))

covMat = Tuple([Symmetric(rand(Float16, i[1], i[2])) for i in zip(nphenotypes, nphenotypes)])
expressionArrays = Tuple([rand() for el in 1:l] for l in ngenes .* ploidy)

parameters = Dict(
  :ngenes => ngenes .* ploidy,
  :nphenotypes => nphenotypes,
  :abiotic_phenotypes => abiotic_phenotypes,
  :growthrates => (0.8, 0.12),
  :pleiotropyMat => pleiotropyMat,
  :epistasisMat =>  epistasisMat,
  :expressionArrays => expressionArrays,
  :selectionCoeffs => (0.5, 0.5),
  :ploidy => ploidy,
  :optPhenotypes => optPhenotypes,
  :covMat => covMat,
  :mutProbs => Tuple([(0.9, 0.0, 0.0), (0.9, 0.0, 0.0)]),
  :mutMagnitudes => Tuple([(0.05, 0.0, 0.01), (0.05, 0.0, 0.01)]),
  :N => Dict(1 => (1000, 1000)),
  :K => Dict(1 => [1000, 1000], 2 => [1000, 1000], 3 => [1000, 1000], 4 => [1000, 1000]),
  # :physical_distance => [[1.0 0.02 0.02 0.02; 0.03 1.0 0.03 0.03; 0.01 0.01 1.0 0.01; 0.01 0.01 0.01 1.0] for i in 1:2],
  :migration_traits => (4, 5),
  :migration_thresholds => (1.5, 1.3),
  :vision_radius => (1,1),
  :check_fraction => (0.5, 0.5),
  :E => (0.01, 0.01),
  :generations => generations,
  :space => space,
  :metric => :chebyshev
)
