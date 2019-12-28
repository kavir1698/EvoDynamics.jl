using Random
import LinearAlgebra: Symmetric


P = (4, 5)  # a tuple  specifying the number of traits p for each species
L = (7, 8)  # a tuple  specifying the number of loci l for each species
m = (2, 1)  # ploidy for each species.
parameters = Dict(
  :L => L .* m,
  :P => P,
  :B =>  Tuple([Random.bitrand(i[1], i[2]) for i in zip(P, L .* m)]),  # a tuple  of pleiotropy matrices, one for each species. Each matrix consists of zeros and ones only. make sure no rows are all zero (a trait is not controled :by any locus)
  :γ => (-0.5, -0.5), # a tuple  of selection coefficients for each species
  :m => m,
  :T => Tuple([randn(Float16, n) for n in P]), # a tuple  of arrays, each inner array θ specifying optimal phenotypes for each species
  :Ω => Tuple([Symmetric(rand(Float16, i[1], i[2])) for i in zip(P, P)]), # a tuple  of matrices, each of which ω represents a covariance matrix of the selection surface
  :M => (0.02, 0.02), # a tuple of mutation rates μ for each species
  :MB => (0.05, 0.05), # a tuple of mutation rates μ of B for each species
  :N => Dict(1 => (1000, 1000)), # a dictionary where a key is node number and the value is a tuple for population size of each species at that node
  :node_capacities => Dict(1 => 2000, 2 => 2000, 3 => 2000, 4 => 2000), # a dictionary where a key is node number and a value is an integer for carrying capacity of the node
  :migration_rates => [1.0 0.02 0.02 0.02; 0.03 1.0 0.03 0.03; 0.01 0.01 1.0 0.01; 0.01 0.01 0.01 1.0],  # a matrix of migration rates between each pair of nodes. The rows and columns of the matrix are node numbers in order.
  :Y => Tuple([rand(Float16, i) for i in L .* m]), # a tuple  of Arrays, each specifying the initial y vector of each species
  :E => (0.8, 0.8), # a tuple  of the variance of a normal distribution ε representing environmental noise for each species.
  :generations => 5, # number of generations to run the simulation
  :space => (2,2),  # Either a tuple for a grid size or a SimpleGraph
  :moore => false
)
