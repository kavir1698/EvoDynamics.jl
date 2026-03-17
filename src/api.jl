"""
    mean_fitness_per_species(model::ABM)

Returns a tuple whose entries are the mean fitness of each species.
"""
function mean_fitness_per_species(model::ABM)
  mean_fitness = Array{Float64}(undef, model.nspecies)
  for species in 1:model.nspecies
    fitness = mean([i.W for i in values(model.agents) if i.species == species])
    mean_fitness[species] = fitness
  end

  return mean_fitness
end

"""
    species_N(model::ABM)

Returns the population size per species in an array.
"""
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
    runmodel(param_file::AbstractString; kwargs)

Creates and runs a model given `parameters`. Returns a `DataFrame` of collected data, which are specified by `kwargs`.

# Keywords

* adata=[] agent data to be collected. Either agent fields or functions that accept an agent as input can be put in the array. To aggregate collected data, provide tuples inside the array. For example, to collect mean and median fitness of individuals which is in field `W`, your array will be [(:W,mean), (:W,median)].
* mdata=[mean_fitness_per_species, species_N] model data to be collected. By default, collects mean population fitness per species. Each row of the output DataFrame corresponds to all agents and each column is the value function applied to a field. The functions in the array are applied to the model object.
* when_data=nothing The generations from which agent data are collected. By default collect at all steps.
* when_model=nothing The generations from which model data are collected. By default collect at all steps.
* replicates::Int = 0 Number of replicates per simulation.
* parallel::Bool = false Whether to run replicates in parallel. If `true`, you should add processors to your julia session (e.g. by `addprocs(n)`) and define your parameters and `EvoDynamics` on all workers. To do that, add `@everywhere` before them. For example, `@everywhere EvoDynamics`.
* seeds = optionally, provide an array of integers as seeds for each replicate.
* agentstep=EvoDynamics.agent_step! Define your own agent stepping if you wish to change the sequence of events or change any one event.
* modelstep=EvoDynamics.model_step!
* showprogress::Bool = false Whether to show a progress meter of the simulations.
* offline_run::Bool = false When false, all the output is saved in memory. When true, It will be written to a file in the intervals specified by `writing_interval`. The files in which it writes to are `adata_filename` and `mdata_filename`. It only works when `parallel = false`.
* writing_interval = 1 See `offline_run`.
* mdata_filename = "mdata.csv" See `offline_run`.
* adata_filename = "adata.csv" See `offline_run`.

"""
function runmodel(param_file::AbstractString;
  adata=nothing, mdata=[mean_fitness_per_species, species_N],
  when_data=nothing,
  when_model=nothing,
  replicates::Int=0,
  parallel::Bool=false,
  seeds=nothing,
  agentstep=EvoDynamics.agent_step!,
  modelstep=EvoDynamics.model_step!,
  showprogress=false,
  offline_run=false,
  writing_interval=1,
  mdata_filename="mdata.csv",
  adata_filename="adata.csv"
)

  if !isnothing(seeds) && replicates > 0
    @assert length(seeds) == replicates "If seeds is not nothing, it has to be an array of ints as long as replicates."
  end

  @info "Loading the parameters..."
  dd = load_parameters(param_file)

  if isnothing(when_data)
    when_data = 0:dd[:generations]
  end
  if isnothing(when_model)
    when_model = 0:dd[:generations]
  end

  # create model
  @info "Running the model..."
  if replicates == 0
    model = model_initiation(dd)
    if offline_run
      offline_run!(model, agentstep, modelstep, model.generations, adata=adata, mdata=mdata, when=when_data, when_model=when_model, agents_first=false, showprogress=showprogress, adata_filename=adata_filename, mdata_filename=mdata_filename,writing_interval=writing_interval)
    else
      # run model and collect data
      agdata, modata = run!(model, agentstep, modelstep, model.generations, adata=adata, mdata=mdata, when=when_data, when_model=when_model, agents_first=false, showprogress=showprogress)
      return agdata, modata, [model]
    end
  else
    models = [model_generator(i, seeds, param_file) for i in 1:replicates]

    agdata, modata, models = ensemblerun!(models, agentstep, modelstep, dd[:generations], adata=adata, mdata=mdata, when=when_data, when_model=when_model, parallel=parallel, agents_first=false, showprogress=showprogress)
    return agdata, modata, models
  end
end

function model_generator(counter, seeds, param_file)
  dd = load_parameters(param_file)
  if isnothing(seeds)
    dd[:seed] = nothing
  else
    dd[:seed] = seeds[counter]
  end
  model = model_initiation(dd)
  return model
end

"""
    expand_iterable_columns(agdata, modata)
  
Expand columns that have arrays (multiple values). Inputs are agent data and model data outputs.
"""
function expand_iterable_columns!(agdata, modata)
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
  return agdata, modata
end