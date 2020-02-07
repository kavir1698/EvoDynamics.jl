using EvoDynamics
using BenchmarkTools
using Random
import LinearAlgebra: Symmetric
Random.seed!(10)


P = (4, 5)
L = (7, 8)
m = (2, 1)
# choose random values for epistasis matrix A, but make sure the diagonal is 
# 1.0, meaning that each locus affects itself 100%.
A = Tuple([Random.rand(-0.5:0.01:0.5, i, i) for i in (L .* m)])
for index in 1:length(A)
  for diag in 1:size(A[index], 1)
    A[index][diag, diag] = 1
  end
end

parameters = Dict(
  :L => L .* m,
  :P => P,
  :R => (0.8, 0.12),
  :C => rand(-0.1:0.01:0.1, 2, 2),
  :B => (rand([true, false], P[1], L[1] * m[1]), rand([true, false], P[2], L[2] * m[2])),
  :A =>  A,
  :Q => Tuple([rand() for el in 1:l] for l in L .* m),
  :Y => (-0.5, -0.5),
  :m => m,
  :T => Tuple([randn(Float16, n) for n in P]),
  :Ω => Tuple([Symmetric(rand(Float16, i[1], i[2])) for i in zip(P, P)]),
  :M => Tuple([(0.02, 0.0, 0.0), (0.02, 0.0, 0.0)]),
  :D => Tuple([(0.05, 0.0, 0.01), (0.05, 0.0, 0.01)]),
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