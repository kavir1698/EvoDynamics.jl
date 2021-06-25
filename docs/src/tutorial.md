# Tutorial

## EvoDynamics.jl's basic usage

First, define your model parameters in a YAML file (here, we call it `parameters.yml`). [Simple Wright-Fisher](@ref) and [Predator prey](@ref) have examples of initiation parameters. See [API](@ref) for a description of each parameter.

We can the use the `runmodel` function to create a model from these parameters and run the simulation.

```@docs
runmodel
```

```@example random
using EvoDynamics
agentdata, modeldata, model = runmodel("parameters.yml")
```

For data collection, see [API](@ref).
