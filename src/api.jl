export runmodel

"Returns a tuple whose entries are the mean fitness of each species."
function mean_fitness_per_species(model::ABM)
  nspecies = length(model.properties[:nphenotypes])
  mean_fitness = Array{Float32}(undef, nspecies)
  for species in 1:nspecies
    fitness = mean([i.W for i in values(model.agents) if i.species == species])
    mean_fitness[species] = fitness
  end

  return (mean_fitness)
end

"""
    runmodel(parameters::Dict; kwargs)

Creates and runs a model given `parameters`. Returns a `DataFrame` of collected data, which are specified by `kwargs`.

# Keywords

* adata=[] agent data to be collected. Either agent fields or functions that accept an agent as input can be put in the array. To aggregate collected data, provide tuples inside the array. For example, to collect mean and median fitness of individuals which is in field `W`, your array will be [(:W,mean), (:W,median)].
* mdata=[mean_fitness_per_species] model data to be collected. By default, collects mean population fitness per species. Each row of the output DataFrame corresponds to all agents and each column is the value function applied to a field. The functions in the array are applied to the model object.
* when::AbstractArray{Int}=1:parameters[:generations] The generations from which data are collected
* replicates::Int = 0 Number of replicates per simulation.
* parallel::Bool = false Whether to run replicates in parallel. If `true`, you should add processors to your julia session (e.g. by `addprocs(n)`) and define your parameters and `EvoDynamics` on all workers. To do that, add `@everywhere` before them. For example, `@everywhere EvoDynamics`.
"""
function runmodel(parameters::Dict;
  adata=nothing, mdata=[mean_fitness_per_species],
  when::AbstractArray{Int}=0:parameters[:generations],
  replicates::Int = 0,
  parallel::Bool = false
  )

  # create model
  model = model_initiation(;parameters...)

  # run model and collect data
  agdata, modata = run!(model, agent_step!, model_step!, parameters[:generations], adata=adata, mdata=mdata, when=when, replicates=replicates, parallel=parallel)

  # Expand columns that have tuples (multiple values)
  allnames = Agents.names(agdata)
  tuple_colsa = Symbol[]
  for nam in allnames
    if eltype(agdata[!, nam]) <: Tuple
      push!(tuple_colsa, nam)
    end
  end
  allnames = Agents.names(modata)
  tuple_colsm = Symbol[]
  for nam in allnames
    if eltype(modata[!, nam]) <: Tuple
      push!(tuple_colsm, nam)
    end
  end

  for nam in tuple_colsa
    nitems = length(agdata[1, nam])
    for item in 1:nitems
      agdata[!, Symbol(string(nam) * "_$item")] = getfield.(agdata[!, nam], item)
    end
    Agents.select!(agdata, Agents.Not(nam))  # remove the column with tuples
  end

  for nam in tuple_colsm
    nitems = length(modata[1, nam])
    for item in 1:nitems
      modata[!, Symbol(string(nam) * "_$item")] = getfield.(modata[!, nam], item)
    end
    Agents.select!(modata, Agents.Not(nam))  # remove the column with tuples
  end

  return agdata, modata, model
end