export runmodel

"Returns a tuple whose entries are the mean fitness of each species."
function mean_fitness_per_species(model::ABM)
  nspecies = length(model.properties[:P])
  mean_fitness = Array{Float32}(undef, nspecies)
  for species in 1:nspecies
    fitness = mean([i.W for i in values(model.agents) if i.species == species])
    mean_fitness[species] = fitness
  end

  return Tuple(mean_fitness)
end

"""
    runmodel(parameters::Dict; kwargs)

Creates and runs a model given `parameters`. Returns a `DataFrame` of collected data, which are specified by `kwargs`.

# Keywords

* collect::Dict=Dict(:model => [mean_fitness_per_species]) Data to be collected. By default, collects mean population fitness per species. Each row of the output DataFrame corresponds to all agents and each column is the value function applied to a field. The functions in a dictionary properties are applied to the collected fields, that is, the keys of properties. For example, to collect mean and median fitness of individuals which is in field `W`, your dictionary will be Dict(:W => [mean, median]).
* when::AbstractArray{Int}=1:parameters[:generations] The generations from which data are collected
* replicates::Int = 0 Number of replicates per simulation.
* parallel::Bool = false Whether to run replicates in parallel. If `true`, you should add processors to your julia session (e.g. by `addprocs(n)`) and define your parameters and `EvoDynamics` on all workers. To do that, add `@everywhere` before them. For example, `@everywhere EvoDynamics`.
"""
function runmodel(parameters::Dict;
  collect::Dict=Dict(:model => [mean_fitness_per_species]),
  when::AbstractArray{Int}=1:parameters[:generations],
  replicates::Int = 0,
  parallel::Bool = false
  )

  # create model
  model = model_initiation(;parameters...)

  # run model and collect data
  data = step!(model, agent_step!, model_step!, parameters[:generations], collect, when=when, replicates=replicates, parallel=parallel)

  # Expand columns that have tuples (multiple values)
  allnames = Agents.names(data)
  tuple_cols = Symbol[]
  for nam in allnames
    if eltype(data[!, nam]) <: Tuple
      push!(tuple_cols, nam)
    end
  end

  for nam in tuple_cols
    nitems = length(data[1, nam])
    for item in 1:nitems
      data[!, Symbol(string(nam) * "_$item")] = getfield.(data[!, nam], item)
    end
    Agents.select!(data, Agents.Not(nam))  # remove the column with tuples
  end

  return data
end