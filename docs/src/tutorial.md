# Tutorial

## EvoDynamics.jl's basic usage

First, define your model parameters. Below is a set of random parameters. See [API](@ref) for a description of each parameter.

```@example random
using Random
import LinearAlgebra: Symmetric

P = (4, 5)
L = (7, 8)
m = (2, 1)
parameters = Dict(
  :L => L .* m,
  :P => P,
  :R => (0.8, 0.12),
  :C => rand(-0.1:0.01:0.1, 2, 2),
  :A =>  Tuple([Random.rand(i[1], i[2]) for i in zip(P, L .* m)]),
  :Y => (-0.5, -0.5),
  :m => m,
  :T => Tuple([randn(Float16, n) for n in P]),
  :Ω => Tuple([Symmetric(rand(Float16, i[1], i[2])) for i in zip(P, P)]),
  :M => (0.02, 0.02),
  :N => Dict(1 => (1000, 1000)),
  :K => Dict(1 => [1000, 1000], 2 => [1000, 1000], 3 => [1000, 1000], 4 => [1000, 1000]),
  :migration_rates => [[1.0 0.02 0.02 0.02; 0.03 1.0 0.03 0.03; 0.01 0.01 1.0 0.01; 0.01 0.01 0.01 1.0] for i in 1:2],
  :E => (0.08, 0.08),
  :generations => 5,
  :space => (2,2),
  :moore => false
)
```

We can the use the `runmodel` function to create a model from these parameters and run the simulation.

```@docs
EvoDynamics.runmodel
```

```@example random
using EvoDynamics
data = runmodel(parameters)
data[1:5, :]
```