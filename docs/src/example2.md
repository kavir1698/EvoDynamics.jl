```@meta
EditURL = "<unknown>/examples/example2.jl"
```

# Weak modular structure

```@example example2
using EvoDynamics
using Agents
using Distributions
using Random
using LinearAlgebra

nspecies = 10
nphenotypes = fill(3, nspecies)
ngenes = fill(3, nspecies)
ploidy = fill(1, nspecies)
epistasisMat = Tuple([rand(Normal(0, 0.01), i, i) for i in (ngenes .* ploidy)])
```

choose random values for epistasis matrix epistasisMat, but make sure the diagonal is
1.0, meaning that each locus affects itself 100%.

```@example example2
for index in 1:length(epistasisMat)
  for diag in 1:size(epistasisMat[index], 1)
    epistasisMat[index][diag, diag] = 1.0
  end
end

mig = rand(10, 10)
mig[LinearAlgebra.diagind(mig)] .= 1.0

parameters = Dict(
  :ngenes => ngenes .* ploidy,
  :nphenotypes => nphenotypes,
  :growthrates => fill(0.1, nspecies),
  :interactionCoeffs => rand(Normal(0, 0.001), nspecies, nspecies),
  :pleiotropyMat => Tuple([LinearAlgebra.diagm(fill(true, max(nphenotypes[i], ngenes[i] * ploidy[i]))) for i in 1:nspecies]),
  :epistasisMat =>  epistasisMat,
  :expressionArrays => Tuple([rand() for el in 1:l] for l in ngenes .* ploidy),
  :selectionCoeffs => fill(0.0, nspecies),
  :ploidy => ploidy,
  :optPhenotypes => Tuple([randn(n) for n in nphenotypes]),
  :covMat => Tuple([LinearAlgebra.diagm(fill(1.0, max(i[1], i[2]))) for i in zip(nphenotypes, nphenotypes)]),
  :mutProbs => Tuple([(0.02, 0.0, 0.0) for i in 1:nspecies]),
  :mutMagnitudes => Tuple([(0.01, 0.0, 0.01) for i in 1:nspecies]),
  :N => Dict(i => fill(100, nspecies) for i in 1:nspecies),
  :K => Dict(i => fill(1000, nspecies) for i in 1:nspecies),
  :migration_rates => [mig for i in 1:nspecies],
  :E => Tuple(0.001 for i in 1:nspecies),
  :generations => 5,
  :space => (5,2),
);

function nspecies_per_node(model)
  output = zeros(model.properties[:nspecies], nv(model))
  for species in model.properties[:nspecies]
    for node in 1:nv(model)
      for ag in ids_in_position(node, model)
        output[model.agents[ag].species, node] += 1
      end
    end
  end
  return Tuple(output)
end

agentdata, modeldata, model = runmodel(parameters, mdata=[nspecies_per_node]);
modeldata
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

