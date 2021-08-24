# Tutorial

## EvoDynamics.jl's basic usage

First, define your model parameters in a YAML file (here, we call it `parameters.yml`). [Simple Wright-Fisher](@ref) and [Predator prey](@ref) have examples of initiation parameters. See [API](@ref) for a description of each parameter.

We can the use the `runmodel` function to create a model from these parameters and run the simulation.

```@docs
runmodel
```

```@example
using EvoDynamics
agentdata, modeldata, model = runmodel("parameters.yml")
```

For data collection, see [API](@ref).

## Creating simulation parameter files

EvoDynamics.jl reads simulation parameters ([API](@ref)) from a human-readable [YAML](https://github.com/JuliaData/YAML.jl) file. This file can be populated manually using a text editor or from within a Julia session.

To create parameters from within a Julia session and write them to a YAML file, you can follow the example below.

To define species parameters, create a dictionary whose keys are numbers starting from 1 for the number of species you have, and values are dictionaries themselves with the species parameter names and values. In the example below, for brevity, I only add a few parameters for each species (see [API](@ref) for the complete list of parameters).


```julia
species_params = Dict(
  1=>Dict("id" => 1, "number of genes" => 2, "number of phenotypes" => 2, "abiotic phenotypes" => [1]),
  2=>Dict("id" => 2, "number of genes" => 2, "number of phenotypes" => 2, "abiotic phenotypes" => [1])
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
