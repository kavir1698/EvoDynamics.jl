# # Host-parasite and predator-prey models

# There are two species interaction equations we can use the Lotka-Voltera competition equation and the generalized Lotka-Voltera equation. The competition equation can model any kind of species interactions except when one species depends on other species to grow, that is, without the other species it has a negative growth rate (e.g. in predator-prey and host-parasite). The generalized equation can model predator-prey and host-parasite dynamics, but it is not suitable for mutualistic dynamics where two species help each other (they may grow indefinitely). 

# Here we model a host-parasite dynamic. The important parameters are `interactionCoeffs`, `interaction_equation`, and `growthrates`. The first species is the parasite that has an intrinsite negative growth rate (-0.1) and completely depends on its host for growth ($ \textrm{interactionCoeffs}_{12} = 0.5 $). The host has a positive intrinsic growth rate (0.1) and it reduces in the presence of parasite. 

using EvoDynamics
using Agents
using Plots

parameters = Dict(
  :ngenes => (2, 2),
  :nphenotypes => (2, 2),
  :growthrates => (-0.1, 0.1),
  :interactionCoeffs => [0 0.5; -0.2 0],
  :interaction_equation => "lotkaVoltera_generalized",
  :pleiotropyMat => ([true false; false true], [true false; false true]),
  :epistasisMat =>  ([1 0; 0 1], [1 0; 0 1]),
  :expressionArrays => ([1, 1], [1, 1]),
  :selectionCoeffs => (0.5, 0.5),
  :ploidy => (1, 2),
  :optPhenotypes => ([2.5, 3], [3.1 2]),
  :covMat => (rand(2, 2), rand(2, 2)),
  :mutProbs => ((0.001, 0.0, 0.0), (0.001, 0.0, 0.0)),
  :mutMagnitudes =>((0.005, 0.0, 0.01), (0.005, 0.0, 0.01)),
  :N => Dict(1 => (100, 100)),
  :K => Dict(1 => (1000, 1000)),
  :migration_rates => [nothing, nothing],
  :E => (0.01, 0.01),
  :generations => 20,
  :space => nothing
)

function nspecies_per_node(model)
  output = zeros(model.properties[:nspecies], nv(model))
  for species in model.properties[:nspecies]
    for node in 1:nv(model)
      for ag in EvoDynamics.ids_in_position(node, model)
        output[model.agents[ag].species, node] += 1
      end
    end
  end
  return Tuple(output)
end

_, modeldata, model = runmodel(parameters, mdata=[nspecies_per_node]);

plot(1:101, modeldata[!, 2], label="Parasite")
plot!(1:101, modeldata[!, 3], label="Host")