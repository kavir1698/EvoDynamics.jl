using EvoDynamics
using BenchmarkTools
using Random
import LinearAlgebra: Symmetric
Random.seed!(10)


nphenotypes = (4, 5)
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

parameters = Dict(
  :ngenes => ngenes .* ploidy,
  :nphenotypes => nphenotypes,
  :growthrates => (0.8, 0.12),
  :interactionCoeffs => rand(-0.1:0.01:0.1, 2, 2),
  :pleiotropyMat => (rand([true, false], nphenotypes[1], ngenes[1] * ploidy[1]), rand([true, false], nphenotypes[2], ngenes[2] * ploidy[2])),
  :epistasisMat =>  epistasisMat,
  :expressionArrays => Tuple([rand() for el in 1:l] for l in ngenes .* ploidy),
  :selectionCoeffs => (0.5, 0.5),
  :ploidy => ploidy,
  :optPhenotypes => Tuple([randn(Float16, n) for n in nphenotypes]),
  :covMat => Tuple([Symmetric(rand(Float16, i[1], i[2])) for i in zip(nphenotypes, nphenotypes)]),
  :mutProbs => Tuple([(0.02, 0.0, 0.0), (0.02, 0.0, 0.0)]),
  :mutMagnitudes => Tuple([(0.05, 0.0, 0.01), (0.05, 0.0, 0.01)]),
  :N => Dict(1 => (1000, 1000)),
  :K => Dict(1 => [1000, 1000], 2 => [1000, 1000], 3 => [1000, 1000], 4 => [1000, 1000]),
  :migration_rates => [[1.0 0.02 0.02 0.02; 0.03 1.0 0.03 0.03; 0.01 0.01 1.0 0.01; 0.01 0.01 0.01 1.0] for i in 1:2],
  :E => (0.01, 0.01),
  :generations => 5,
  :space => (2,2),
  :moore => false
)

a = @benchmark runmodel(parameters)
display(a)