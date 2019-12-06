using Random
import LinearAlgebra: Symmetric

P = (4, 5)  # a tuple  specifying the number of traits p for each species
L = (7, 8)  # a tuple  specifying the number of loci l for each species
m = 1  # ploidy, m=2 diploid. Currently only haploids are implemented
parameters = Dict(
  :L => L,
  :P => P,
  :B =>  Tuple([Random.bitrand(i[1], i[2]) for i in zip(P, L)]),  # a tuple  of pleiotropy matrices, one for each species. Each matrix consists of zeros and ones only. make sure no rows are all zero (a trait is not controled :by any locus)
  :γ => (-0.5, -0.5), # a tuple  of selection coefficients for each species
  :m => m,
  :T => Tuple([randn(Float16, n) for n in P]), # a tuple  of arrays, each inner array θ specifying optimal phenotypes for each species
  :Ω => Tuple([Symmetric(rand(Float16, i[1], i[2])) for i in zip(P, P)]), # a tuple  of matrices, each of wich ω represents a covariance matrix of the selection surface
  :M => (0.02, 0.02), # a tuple of mutation rates μ for each species
  :MB => (0.05, 0.05), # a tuple of mutation rates μ<sub>B</sub> for each species
  :N => (1000, 1000), # a tuple  for population size of each species
  :Y => Tuple([rand(Float16, i*m) for i in L]), # a tuple  of Arrays, each specifying the initial y vector of each species
  :E => (0.8, 0.8), # a tuple  of the variance of a normal distribution ε representing environmental noise for each species.
  :generations => 15 # number of generations to run the simulation
)
