# Tutorial

## EvoDynamics.jl's basic usage

First, define your model parameters in a YAML file (here, we call it `parameters.yml`). [Simple Wright-Fisher](@ref) and [Predator prey](@ref) have examples of initiation parameters. See [Model description](@ref) for a description of each parameter.

We can the use the `runmodel` function to create a model from these parameters and run the simulation.

```@docs
runmodel
```

```@example
using EvoDynamics
agentdata, modeldata, model = runmodel("parameters.yml")
```

For data collection, see [Model description](@ref).

## Creating simulation parameter files

EvoDynamics.jl reads simulation parameters ([Model description](@ref)) from a human-readable [YAML](https://github.com/JuliaData/YAML.jl) file. This file can be populated manually using a text editor or from within a Julia session.

To create parameters from within a Julia session and write them to a YAML file, you can follow the example below.

To define species parameters, create a dictionary whose keys are numbers starting from 1 for the number of species you have, and values are dictionaries themselves with the species parameter names and values. In the example below, for brevity, I only add a few parameters for each species (see [Model description](@ref) for the complete list of parameters).


```julia
species_params = Dict(
  1=>Dict("name" => "a", "number of genes" => 2, "number of phenotypes" => 2, "abiotic phenotypes" => [1]),
  2=>Dict("name" => "b", "number of genes" => 2, "number of phenotypes" => 2, "abiotic phenotypes" => [1])
)

```

To define model parameters, create a new dictionary whose keys are the parameter names and values the parameter values. Here, again I only define a few parameters for brevity.

```julia

model_params = Dict("generations"=> 100, "space"=> [6,10], "food sources" => [1.0, 0.7, 0.0, 0.0])

```

Finally, create a dictionary mixing the two dictionaries before with keys "species" and "model". 

```julia
data = Dict("species"=> species_params, "model"=> model_params)
```

This dictionary can be written as a YAML file in the correct format.

```julia
using YAML
f = "params.yml"
YAML.write_file(f, data)
```


## Data collection

The interface to the model is from the `runmodel` function (see [Tutorial](@ref)).

EvoDynamics.jl uses [Agents.jl](https://github.com/JuliaDynamics/Agents.jl) underneath. See [Agents.jl's documentation](https://juliadynamics.github.io/Agents.jl/dev/) for details about writing functions to collect any data during simulations. Here, we explain the specific implementation of the model.

There are two main objects from which you can collect data: and agent object of type `AbstractAgent` and a model object of type `ABM`. Both of these types are defined the `Agents.jl` package.

Agent object has the following fields that are individual specific: `id`, `pos`, `species`, `epistasisMat` (epistasis matrix), `pleiotropyMat` (pleiotropy matrix), `q` (gene expression array), `biotic_phenotype`, `abiotic_phenotype`, `age`, `sex`, `energy`, `interaction_history`, and `W` (fitness).

The model object has the following fields that can be accessed with the `.` syntax and describe properties of species or the model: `ngenes`, `nphenotypes`, `growthrates`, `selectionCoeffs`, `ploidy`, `optvals` (optimal values), `optinds` (optval indices per generation), `mutProbs` (mutation probabilities), `mutMagnitudes` (mutation magnitudes), `N`, `E` (environmental noise), `generations`, `nspecies`, `migration_traits`, `vision_radius`, `check_fraction`, `migration_thresholds`, `step`, `biotic_phenotypes` (indices of biotic phenotypes per species), `abiotic_phenotypes`, `max_ages`, `food_sources`, `interactions`, `resources`, `recombination`, `initial_energy`.

You can collect data from agents and/or from the model object. To collect data from agents, use the `adata` keyword argument in the `runmodel` function, and to collect data from the model, use the `mdata` keyword. A complete description of the values these keywords take are at [data collection section of the Agents.jl package](https://juliadynamics.github.io/Agents.jl/stable/tutorial/#.-Collecting-data).

For example, we use the function below to count the number of individual per species:

```jl
"Returns the size of each species."
function species_N(model::ABM)
  allagents = model.agents
  if length(allagents) == 0
    return fill(0, model.nspecies)
  else
    counts = countmap([a.species for a in values(model.agents)])
    output = fill(0, model.nspecies)
    for (k, v) in counts
      output[k] = v
    end
    return output
  end
end

using EvoDynamics

agentdata, modeldata, model = runmodel("parameters.yml", mdata=[species_N])
```
