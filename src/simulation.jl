
"""
A `struct` for individuals that keeps individual-specific variables.
"""
mutable struct Ind{B<:AbstractFloat, C<:AbstractArray, D<:AbstractArray, E<:AbstractArray} <: AbstractAgent
  id::Int  # the individual ID
  pos::Tuple{Int, Int}  # the individuals position
  species::Int  # the species ID the individual belongs to
  biotic_phenotype::Vector{B}
  abiotic_phenotype::Vector{B}
  epistasisMat::C  # epistasis matrix
  pleiotropyMat::D  # pleiotropy matrix
  q::E  # expression array
  age::Int
  sex::Bool
  interaction_history::MArray{S, Int} where S # records the last interaction with all species
  energy::B  # determines whether need to feed (energy=0) or not.
  W::B  # survival probability
end

struct Params{F<:AbstractFloat, I<:Int}
  ngenes::Vector{I}
  nphenotypes::Vector{I}
  growthrates::Vector{F}
  selectionCoeffs::Vector{F}
  ploidy::Vector{I}
  optvals::Vector{Vector{Vector{Matrix{F}}}}
  optinds::Vector{Vector{I}}
  mutProbs::Vector{Vector{DiscreteNonParametric{Bool, Float64, Vector{Bool}, Vector{Float64}}}}
  mutMagnitudes::Vector{Vector{UnivariateDistribution{S} where S<:ValueSupport}}
  N::Vector{Vector{I}}
  E::Vector{Normal{F}}
  generations::I
  nspecies::I
  new_N::Vector{Vector{I}}
  migration_traits::Vector{I}
  vision_radius::Vector{I}
  check_fraction::Vector{F}
  migration_thresholds::Vector{F}
  step::MVector{1, Int64}
  nodes::Matrix{Tuple{I, I}}
  biotic_phenotypes::Vector{Vector{I}}
  abiotic_phenotypes::Vector{Vector{I}}
  max_ages::Vector{I}
  ids::Dict{I, I}
  food_sources::Matrix{F}
  interactions::Matrix{F}
  resources::Matrix{I}
  resources_org::Matrix{I}
  recombination::Vector{Poisson{F}}
  initial_energy::Vector{F}
end

const variance = 1.0

"""
    model_initiation(param_file)

Innitializes the model.
"""
function model_initiation(param_file)
  dd = load_parameters(param_file)

  if !isnothing(dd[:model]["seed"])
    Random.seed!(dd[:model]["seed"])
  end
  
  properties, species_arrays = create_properties(dd)
  epistasisMat, pleiotropyMat, expressionArrays = species_arrays

  space = dd[:model]["space"]

  if isnothing(space)
    fspace = GridSpace((1, 1))
  elseif typeof(space) <: Tuple
    fspace = GridSpace(space, periodic=dd[:model]["periodic"], metric=dd[:model]["metric"])
  end

  indtype = EvoDynamics.Ind{typeof(0.1), eltype(epistasisMat), eltype(pleiotropyMat), eltype(expressionArrays)}

  model = ABM(indtype, fspace, properties=properties, scheduler=Schedulers.randomly)
  
  # create and add agents
  for sp in 1:model.nspecies
    for (pos, n) in enumerate(model.N[sp])
      for ind in 1:n
        abiotic_ph = get_abiotic_phenotype(sp, epistasisMat[sp], pleiotropyMat[sp], expressionArrays[sp], model)
        biotic_ph = get_biotic_phenotype(sp, epistasisMat[sp], pleiotropyMat[sp], expressionArrays[sp], model)
        W = abiotic_fitness(abiotic_ph, sp, pos, model)
        sex = false
        if model.ploidy[sp] == 2
          sex = rand((true, false))
        end
        interaction_history = MVector{model.nspecies, Int}(fill(0, model.nspecies))
        initial_energy = model.initial_energy[sp]
        add_agent!(model.nodes[pos], model, sp, biotic_ph, abiotic_ph, epistasisMat[sp], pleiotropyMat[sp], expressionArrays[sp], 0, sex, interaction_history, initial_energy, W)
      end      
    end
  end

  return model
end

nnodes(x::GridSpace) = prod(size(x))
nnodes(x::ABM) = nnodes(x.space)

function create_properties(dd)
  nspecies = length(dd[:species])

  Ed = [Normal(0.0, dd[:species][i]["environmental noise"]) for i in 1:nspecies]
  Mdists = [[DiscreteNonParametric([true, false], [i, 1-i]) for i in dd[:species][arr]["mutation probabilities"]] for arr in 1:nspecies]  # μ (probability of change)
  Ddists = [[Normal(0, dd[:species][ar]["mutation magnitudes"][1]), DiscreteNonParametric([true, false], [dd[:species][ar]["mutation magnitudes"][2], 1-dd[:species][ar]["mutation magnitudes"][2]]), Normal(0, dd[:species][ar]["mutation magnitudes"][3])] for ar in 1:nspecies]  # amount of change in case of mutation
  
  # make single-element arrays 2D so that linAlg functions will work
  newA = Array{Array{Float64}}(undef, nspecies)
  newQ = Array{Array{Float64}}(undef, nspecies)
  for i in 1:nspecies
    if length(dd[:species][i]["epistasis matrix"]) == 1
      newA[i] = reshape(dd[:species][i]["epistasis matrix"], 1, 1)
      newQ[i] = reshape(dd[:species][i]["expression array"], 1, 1)
    else
      newA[i] = dd[:species][i]["epistasis matrix"]
      newQ[i] = dd[:species][i]["expression array"]
    end
  end

  epistasisMatS = [MArray{Tuple{size(newA[i])...}}(newA[i]) for i in eachindex(newA)]
  pleiotropyMatS = [MArray{Tuple{size(dd[:species][i]["pleiotropy matrix"])...}}(dd[:species][i]["pleiotropy matrix"]) for i in 1:nspecies]
  expressionArraysS = [MArray{Tuple{size(newQ[i])...}}(newQ[i]) for i in eachindex(newQ)]

  ngenes = [dd[:species][i]["number of genes"] for i in 1:nspecies]
  nphenotypes = [dd[:species][i]["number of phenotypes"] for i in 1:nspecies]
  growthrates = [dd[:species][i]["growth rate"] for i in 1:nspecies]
  selectionCoeffs = [dd[:species][i]["selection coefficient"] for i in 1:nspecies]
  ploidy = [dd[:species][i]["ploidy"] for i in 1:nspecies]
  optvals = [dd[:species][i]["optimal phenotype values"] for i in 1:nspecies]
  optinds = [dd[:species][i]["optimal phenotypes"] for i in 1:nspecies]
  Ns = [dd[:species][i]["N"] for i in 1:nspecies]
  migration_traits = [dd[:species][i]["migration phenotype"] for i in 1:nspecies]
  vision_radius = [dd[:species][i]["vision radius"] for i in 1:nspecies]
  check_fraction = [dd[:species][i]["check fraction"] for i in 1:nspecies]
  migration_thresholds = [dd[:species][i]["migration threshold"] for i in 1:nspecies]
  generations = dd[:model]["generations"]
  step = MVector{1, Int}(undef)
  step[1] = 0
  nnodes = Matrix{Tuple{Int, Int}}(undef, dd[:model]["space"]...)
  for col in 1:size(nnodes, 2)
    for row in 1:size(nnodes, 1)
      nnodes[row, col] = (row, col)
    end
  end
  biotic_phenotyps = [dd[:species][i]["biotic phenotypes"] for i in 1:nspecies]
  abiotic_phenotyps = [dd[:species][i]["abiotic phenotypes"] for i in 1:nspecies]
  max_ages = [dd[:species][i]["age"] for i in 1:nspecies]
  ids = Dict(i => dd[:species][i]["id"] for i in 1:nspecies)
  recombination = [Poisson(dd[:species][i]["recombination"]) for i in 1:nspecies]
  initial_energy = [AbstractFloat(dd[:species][i]["initial energy"]) for i in 1:nspecies]

  properties = Params(ngenes, nphenotypes, growthrates, selectionCoeffs, ploidy, optvals, optinds, Mdists, Ddists, Ns, Ed, generations, nspecies, Ns, migration_traits, vision_radius, check_fraction, migration_thresholds, step, nnodes, biotic_phenotyps, abiotic_phenotyps, max_ages, ids, dd[:model]["food sources"], dd[:model]["interactions"], dd[:model]["resources"], deepcopy(dd[:model]["resources"]), recombination, initial_energy)
  
  return properties, (epistasisMatS, pleiotropyMatS, expressionArraysS)
end

function return_opt_phenotype(species::Int, generation::Int, site::Int, model::ABM)
  nabiotic = length(model.abiotic_phenotypes[species])
  output = Array{typeof(0.1)}(undef, nabiotic)
  for trait in 1:nabiotic
    output[trait] = model.optvals[species][model.optinds[species][generation+1]][trait][site]
  end
  return output
end

function return_opt_phenotype(species::Int, generation::Int, site::Tuple{Int, Int}, model::ABM)
  nabiotic = length(model.abiotic_phenotypes[species])
  output = Array{typeof(0.1)}(undef, nabiotic)
  for trait in 1:nabiotic
    output[trait] = model.optvals[species][model.optinds[species][generation+1]][trait][site[1],site[2]]
  end
  return output
end

function model_step!(model::ABM)
  model.step[1] += 1
  model.resources .= model.resources_org
end

function agent_step!(agent::Ind, model::ABM)
  # update age
  agent.age += 1
  # migrate
  migrate!(agent, model)
  # use food
  use_energy!(agent, model)
  # consume basic energy if agent can
  consume_food!(agent, model)
  # interact with other species
  interact!(agent, model)
  # reproduction for the haploid
  reproduce!(agent, model)
  # survive
  survive!(agent, model)
end

"""
    survive!(agent::Ind, model::ABM)

Kills the agent if it has no energy, is too old, or by change given its fitness `W`
"""
function survive!(agent::Ind, model::ABM)
  if !haskey(model.agents, agent.id)
    return
  elseif agent.energy < 0
    kill_agent!(agent, model)
  elseif agent.age ≥ model.max_ages[agent.species]
    kill_agent!(agent, model)
  elseif rand() > adjusted_fitness(agent, model)
    kill_agent!(agent, model)
  end
end

function adjusted_fitness(agent, model)
  W = agent.W < 0 ? 0.0 : agent.W
  1.0 - ( (1.0 - W) * model.selectionCoeffs[agent.species])
end

function use_energy!(agent, model)
  agent.energy -= sum(model.food_sources[agent.species, :])
end

function consume_food!(agent::Ind, model::ABM)
  if model.food_sources[agent.species, agent.species] > 0 && model.resources[agent.pos...] > 0
    model.resources[agent.pos...] -= 1
    agent.energy += 1
  end
end

"Mutate an agent."
function mutate!(agent::Ind, model::ABM)
  mutated = false
  # mutate gene expression
  if rand(model.mutProbs[agent.species][1])
    agent.q .+= rand(model.mutMagnitudes[agent.species][1], model.ngenes[agent.species] * model.ploidy[agent.species])
    mutated = true
  end
  # mutate pleiotropy matrix
  if rand(model.mutProbs[agent.species][2])
    randnumbers = rand(model.mutMagnitudes[agent.species][2], size(agent.pleiotropyMat))
    agent.pleiotropyMat[randnumbers] .= .!agent.pleiotropyMat[randnumbers]
    mutated = true
  end
  # mutate epistasis matrix
  if rand(model.mutProbs[agent.species][3])
    agent.epistasisMat .+= rand(model.mutMagnitudes[agent.species][3], size(agent.epistasisMat))
    mutated = true
  end
  # update biotic and abiotic phenotypes and W
  if mutated
    abiotic_phenotype = get_abiotic_phenotype(agent.species, agent.epistasisMat, agent.pleiotropyMat, agent.q, model) 
    biotic_phenotype = get_biotic_phenotype(agent.species, agent.epistasisMat, agent.pleiotropyMat, agent.q, model)
    W = abiotic_fitness(abiotic_phenotype, agent.species, agent.pos, model)
    agent.biotic_phenotype .= biotic_phenotype
    agent.abiotic_phenotype .= abiotic_phenotype
    agent.W = W
  end
end

# coord2vertex(c, model) = c[1] + (size(model.space)[1] * (c[2]-1))
