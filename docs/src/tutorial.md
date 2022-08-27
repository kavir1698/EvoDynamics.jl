# Tutorial

## Installation

Install using the following command in a Julia REPL.

```julia
]add EvoDynamics
```

## Basic usage

Parameters of a model should be put in a julia file (.jl format) with the structure below.

See [Simple Wright-Fisher](@ref) and [Predator prey](@ref) for complete examples of parameter files.

```jl
## 1. Functions
...

## 2. Species parameters. A dictionary for each species

## 3. Model parameters as a dictionary.

```

The order of these sections is important because each section uses objects from its preceding sections. 

Functions are used to create parameters that may change temporally and spatially. The following parameters are functions: bottleneck function which kills certain agents at certain times and spaces, optimal phenotype values which returns the optimal phenotype for a species at a given time and space, and environmental resources that may change over time.


You may create as many species as you want. Parameters of each species is a dictionary.

Model parameters is one dictionary that stores general parameters of the model, such as number of generations, space size, and species interaction parameters.

First, define your model parameters (here, we call it `parameters.jl`). [Simple Wright-Fisher](@ref) and [Predator prey](@ref) have examples of initiation parameters. See [Model description](@ref) for a description of each parameter.

We can the use the `runmodel` function to create a model from these parameters and run the simulation.

```@docs
runmodel
```

```@example
using EvoDynamics
agentdata, modeldata, model = runmodel("parameters.jl")
```

## Creating simulation parameter files

EvoDynamics.jl reads simulation parameters ([Model description](@ref)) from a julia file containing dictionaries and functions. This file can be populated manually using a text editor or from within a Julia session.

## Collecting data

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

## Running simulations in parallel

You can run replicate simulation in parallel. To that end, you need to add processors, and import EvoDynamics.jl and your parameters files on all cores:

```julia
using Distributed
addprocs(4)
@everywhere using EvoDynamics
@everywhere param_file = "params.jl"
@everywhere EvoDynamics.load_parameters(param_file)
adata, mdata, models = runmodel(param_file, replicates=10, parallel=true)
```