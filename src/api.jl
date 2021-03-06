export runmodel

"Returns a tuple whose entries are the mean fitness of each species."
function mean_fitness_per_species(model::ABM)
  mean_fitness = Array{Float64}(undef, model.nspecies)
  for species in 1:model.nspecies
    fitness = mean([i.W for i in values(model.agents) if i.species == species])
    mean_fitness[species] = fitness
  end

  return mean_fitness
end

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

"""
    runmodel(parameters::Dict; kwargs)

Creates and runs a model given `parameters`. Returns a `DataFrame` of collected data, which are specified by `kwargs`.

# Keywords

* adata=[] agent data to be collected. Either agent fields or functions that accept an agent as input can be put in the array. To aggregate collected data, provide tuples inside the array. For example, to collect mean and median fitness of individuals which is in field `W`, your array will be [(:W,mean), (:W,median)].
* mdata=[mean_fitness_per_species] model data to be collected. By default, collects mean population fitness per species. Each row of the output DataFrame corresponds to all agents and each column is the value function applied to a field. The functions in the array are applied to the model object.
* when=nothing The generations from which data are collected. By default collect at all steps.
* replicates::Int = 0 Number of replicates per simulation.
* parallel::Bool = false Whether to run replicates in parallel. If `true`, you should add processors to your julia session (e.g. by `addprocs(n)`) and define your parameters and `EvoDynamics` on all workers. To do that, add `@everywhere` before them. For example, `@everywhere EvoDynamics`.
"""
function runmodel(param_file;
  adata=nothing, mdata=[mean_fitness_per_species, species_N],
  when=nothing,
  replicates::Int = 0,
  parallel::Bool = false
  )

  # create model
  model = model_initiation(param_file)

  if isnothing(when)
    whenn = 0:model.generations
  else
    whenn = when
  end

  # run model and collect data
  agdata, modata = run!(model, agent_step!, model_step!, model.generations, adata=adata, mdata=mdata, when=whenn, replicates=replicates, parallel=parallel, agents_first=false)

  # Expand columns that have arrays (multiple values)
  allnames = Agents.names(agdata)
  tuple_colsa = String[]
  for nam in allnames
    if eltype(agdata[!, nam]) <: AbstractArray
      push!(tuple_colsa, nam)
    end
  end
  allnames = Agents.names(modata)
  tuple_colsm = String[]
  for nam in allnames
    if eltype(modata[!, nam]) <: AbstractArray
      push!(tuple_colsm, nam)
    end
  end

  for nam in tuple_colsa
    nitems = length(agdata[1, nam])
    for item in 1:nitems
      agdata[!, Symbol(nam * "_$item")] = getindex.(agdata[!, Symbol(nam)], item)
    end
    Agents.select!(agdata, Agents.Not(nam))  # remove the column with tuples
  end

  for nam in tuple_colsm
    nitems = length(modata[1, nam])
    for item in 1:nitems
      modata[!, Symbol(nam * "_$item")] = getindex.(modata[!, Symbol(nam)], item)
    end
    Agents.select!(modata, Agents.Not(Symbol(nam)))  # remove the column with tuples
  end

  return agdata, modata, model
end