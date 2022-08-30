
"""
A `struct` for individuals that keeps individual-specific variables.
"""
mutable struct Ind{B<:AbstractFloat,C<:AbstractArray,D<:AbstractArray,E<:AbstractArray} <: AbstractAgent
  id::Int  # the individual ID
  pos::Tuple{Int,Int}  # the individuals position
  species::Int  # the species ID the individual belongs to
  biotic_phenotype::Vector{B}
  abiotic_phenotype::Vector{B}
  epistasisMat::C  # epistasis matrix
  pleiotropyMat::D  # pleiotropy matrix
  q::E  # expression array
  age::Int
  sex::Bool
  interaction_history::MArray{S,Int} where {S} # records the last interaction with all species
  energy::B  # determines whether need to feed (energy=0) or not.
  W::B  # survival probability
  isalive::Bool
  mate::Int  # id of mate. Only relevant for diploids
  time_met_other_sex::Int
end

struct Params{F<:AbstractFloat,I<:Int,N<:AbstractString}
  ngenes::Vector{I}
  nphenotypes::Vector{I}
  growthrates::Vector{F}
  selectionCoeffs::Vector{F}
  ploidy::Vector{I}
  optvals::Vector{Vector{Matrix{Vector{F}}}}
  mutProbs::Vector{Vector{DiscreteNonParametric{Bool,Float64,Vector{Bool},Vector{Float64}}}}
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
  step::MVector{1,Int64}
  nodes::Matrix{Tuple{I,I}}
  biotic_phenotypes::Vector{Vector{I}}
  abiotic_phenotypes::Vector{Vector{I}}
  max_ages::Vector{I}
  names::Dict{I,N}
  food_sources::Matrix{F}
  interactions::Matrix{F}
  resources::Matrix{I}
  resources_org::Vector{Matrix{I}}
  recombination::Vector{Poisson{F}}
  initial_energy::Vector{F}
  bottlenecks::Vector{Vector{Matrix{F}}}
  repro_start::Vector{Int}
  repro_end::Vector{Int}
  biotic_variances::Vector{F}
  abiotic_variances::Vector{F}
  mating_schemes::Vector{Int}
end

"""
    model_initiation(dd)

Innitializes the model.
"""
function model_initiation(dd)

  if isnothing(dd[:seed])
    mrng = RandomNumbers.Xorshifts.Xoroshiro128Plus()
  else
    mrng = RandomNumbers.Xorshifts.Xoroshiro128Plus(dd[:seed])
  end

  if isnothing(dd[:space])
    dd[:space] = (1, 1)
  end

  properties, species_arrays = create_properties(dd)
  epistasisMat, pleiotropyMat, expressionArrays = species_arrays

  space = dd[:space]

  if isnothing(space)
    fspace = GridSpace((1, 1))
  elseif typeof(space) <: Tuple
    fspace = GridSpace(space, periodic=dd[:periodic], metric=Symbol(dd[:metric]))
  end

  indtype = EvoDynamics.Ind{typeof(0.1),eltype(epistasisMat),eltype(pleiotropyMat),eltype(expressionArrays)}

  model = ABM(indtype, fspace, properties=properties, scheduler=Schedulers.randomly, rng=mrng)

  # create and add agents
  for sp in 1:model.nspecies
    for (pos, n) in enumerate(model.N[sp])
      for ind in 1:n
        abiotic_ph = get_abiotic_phenotype(sp, epistasisMat[sp], pleiotropyMat[sp], expressionArrays[sp], model)
        biotic_ph = get_biotic_phenotype(sp, epistasisMat[sp], pleiotropyMat[sp], expressionArrays[sp], model)
        W = abiotic_fitness(abiotic_ph, sp, pos, model)
        W = adjust_fitness(W, sp, model)
        sex = false
        if model.ploidy[sp] == 2
          sex = rand(model.rng, (true, false))
        end
        interaction_history = MVector{model.nspecies,Int}(fill(-1, model.nspecies))
        initial_energy = model.initial_energy[sp]
        add_agent!(model.nodes[pos], model, sp, biotic_ph, abiotic_ph, epistasisMat[sp], pleiotropyMat[sp], expressionArrays[sp], 1, sex, interaction_history, initial_energy, W, true, 0, -1)
      end
    end
  end

  return model
end

nnodes(x::GridSpace) = prod(size(x))
nnodes(x::ABM) = nnodes(x.space)

function create_properties(dd)
  nspecies = length(dd[:species])
  allspecies = dd[:species]
  Ed = [Normal(0.0, allspecies[i][:environmental_noise]) for i in 1:nspecies]
  Mdists = [[DiscreteNonParametric([true, false], [i, 1 - i]) for i in allspecies[arr][:mutation_probabilities]] for arr in 1:nspecies]  # μ (probability of change)
  Ddists = [[Normal(0, allspecies[ar][:mutation_magnitudes][1]), DiscreteNonParametric([true, false], [allspecies[ar][:mutation_magnitudes][2], 1 - allspecies[ar][:mutation_magnitudes][2]]), Normal(0, allspecies[ar][:mutation_magnitudes][3])] for ar in 1:nspecies]  # amount of change in case of mutation

  # make single-element arrays 2D so that linAlg functions will work
  newA = Array{Array{Float64}}(undef, nspecies)
  newQ = Array{Array{Float64}}(undef, nspecies)
  for i in 1:nspecies
    if length(allspecies[i][:epistasis_matrix]) == 1
      newA[i] = reshape(allspecies[i][:epistasis_matrix], 1, 1)
      newQ[i] = reshape(allspecies[i][:expression_array], 1, 1)
    else
      newA[i] = allspecies[i][:epistasis_matrix]
      newQ[i] = allspecies[i][:expression_array]
    end
  end

  epistasisMatS = [MArray{Tuple{size(newA[i])...}}(newA[i]) for i in eachindex(newA)]
  pleiotropyMatS = [MArray{Tuple{size(allspecies[i][:pleiotropy_matrix])...}}(Bool.(allspecies[i][:pleiotropy_matrix])) for i in 1:nspecies]
  expressionArraysS = [MArray{Tuple{size(newQ[i])...}}(newQ[i]) for i in eachindex(newQ)]

  ngenes = [allspecies[i][:number_of_genes] for i in 1:nspecies]
  nphenotypes = [allspecies[i][:number_of_phenotypes] for i in 1:nspecies]
  growthrates = [allspecies[i][:growth_rate] for i in 1:nspecies]
  selectionCoeffs = [allspecies[i][:selection_coefficient] for i in 1:nspecies]
  ploidy = [allspecies[i][:ploidy] for i in 1:nspecies]
  optvals = [allspecies[i][:optimal_phenotypes] for i in 1:nspecies]
  Ns = [allspecies[i][:N] for i in 1:nspecies]
  migration_traits = [allspecies[i][:migration_phenotype] for i in 1:nspecies]
  vision_radius = [allspecies[i][:vision_radius] for i in 1:nspecies]
  check_fraction = [allspecies[i][:check_fraction] for i in 1:nspecies]
  migration_thresholds = [allspecies[i][:migration_threshold] for i in 1:nspecies]
  generations = dd[:generations]
  step = MVector{1,Int}(undef)
  step[1] = 0
  nnodes = Matrix{Tuple{Int,Int}}(undef, dd[:space]...)
  for col in 1:size(nnodes, 2)
    for row in 1:size(nnodes, 1)
      nnodes[row, col] = (row, col)
    end
  end
  biotic_phenotyps = [allspecies[i][:biotic_phenotypes] for i in 1:nspecies]
  abiotic_phenotyps = [allspecies[i][:abiotic_phenotypes] for i in 1:nspecies]
  max_ages = [allspecies[i][:age] for i in 1:nspecies]
  names = Dict(i => allspecies[i][:name] for i in 1:nspecies)
  recombination = [Poisson(allspecies[i][:recombination]) for i in 1:nspecies]
  initial_energy = [AbstractFloat(allspecies[i][:initial_energy]) for i in 1:nspecies]
  bottlenecks = [allspecies[i][:bottleneck_function] for i in 1:nspecies]
  repro_start = [allspecies[i][:reproduction_start_age] for i in 1:nspecies]
  repro_end = [allspecies[i][:reproduction_end_age] for i in 1:nspecies]
  resources = dd[:resources][1]
  biotic_variances = [allspecies[i][:biotic_variance] for i in 1:nspecies]
  abiotic_variances = [allspecies[i][:abiotic_variance] for i in 1:nspecies]
  mating_schemes = [allspecies[i][:mating_scheme] for i in 1:nspecies]

  # reshape single value matrices to (1,1)
  if length(resources) == 1
    resources = reshape(resources, 1, 1)
    dd[:food_sources] = reshape(dd[:food_sources], 1, 1)
    dd[:interactions] = reshape(dd[:interactions], 1, 1)
  end

  properties = Params(ngenes, nphenotypes, growthrates, selectionCoeffs, ploidy, optvals, Mdists, Ddists, Ns, Ed, generations, nspecies, Ns, migration_traits, vision_radius, check_fraction, migration_thresholds, step, nnodes, biotic_phenotyps, abiotic_phenotyps, max_ages, names, dd[:food_sources], dd[:interactions], resources, dd[:resources], recombination, initial_energy, bottlenecks, repro_start, repro_end, biotic_variances, abiotic_variances, mating_schemes)

  return properties, (epistasisMatS, pleiotropyMatS, expressionArraysS)
end

function return_opt_phenotype(species::Int, site::Int, model::ABM)
  # ss = vertex2coord(site, model)
  model.optvals[species][model.step[1]+1][site]
end

function return_opt_phenotype(species::Int, site::Tuple{Int,Int}, model::ABM)
  model.optvals[species][model.step[1]+1][site...]
end

function model_step!(model::ABM)
  model.step[1] += 1
  model.resources .= model.resources_org[model.step[1]+1]
end

function agent_step!(agent::Ind, model::ABM)
  # abiotic survive
  abiotic_survive!(agent, model)
  if !agent.isalive
    return
  end
  # use food
  burn_energy!(agent)
  # consume basic energy if agent can
  consume_food!(agent, model)
  # interact with other species: predation, cooperation, competition.
  interact!(agent, model)
  # Kill the agent if it doesn't have energy
  if agent.isalive && agent.energy < 0
    remove_agent!(agent, model)
    return
  end
  if !agent.isalive
    return
  end
  # reproduction for both haploids and diploids
  reproduce!(agent, model)
  # migrate
  migrate!(agent, model)

  if agent.isalive && agent.age ≥ model.max_ages[agent.species]
    remove_agent!(agent, model)
  end
  # bottleneck
  bn_prob = model.bottlenecks[agent.species][model.step[1]+1][agent.pos...]
  if agent.isalive && bn_prob > 0.0
    if rand(model.rng) < bn_prob
      remove_agent!(agent, model)
    end
  end
  # update age
  agent.age += 1
end

function remove_agent!(agent, model)
  agent.isalive = false
  kill_agent!(agent, model)
end

"""
    survive!(agent::Ind, model::ABM)

Kills the agent if it has no energy, or is too old.
"""
function survive!(agent::Ind, model::ABM)
  if !agent.isalive
    return
  elseif agent.energy < 0
    remove_agent!(agent, model)
  elseif agent.age ≥ model.max_ages[agent.species]
    remove_agent!(agent, model)
  end
end

"""
Kills the agent by chance given its fitness `W`.
"""
function abiotic_survive!(agent::Ind, model::ABM)
  if agent.isalive && rand(model.rng) > agent.W
    remove_agent!(agent, model)
  end
end

function adjust_fitness!(agent::Ind, model::ABM)
  W = agent.W < 0 ? 0.0 : agent.W
  newW = 1.0 - ((1.0 - W) * model.selectionCoeffs[agent.species])
  agent.W = newW
end

function adjust_fitness(W, species, model::ABM)
  W2 = W < 0 ? 0.0 : W
  newW = 1.0 - ((1.0 - W2) * model.selectionCoeffs[species])
  return newW
end

function burn_energy!(agent)
  agent.energy -= 1
end

function consume_food!(agent::Ind, model::ABM)
  environmental_consumption = model.food_sources[agent.species, agent.species]
  if environmental_consumption > 0 && model.resources[agent.pos...] >= environmental_consumption
    model.resources[agent.pos...] -= environmental_consumption
    agent.energy += environmental_consumption
  end
end

"Mutate an agent."
function mutate!(agent::Ind, model::ABM)
  mutated = false
  # mutate gene expression
  if rand(model.rng, model.mutProbs[agent.species][1])
    agent.q .+= rand(model.rng, model.mutMagnitudes[agent.species][1], model.ngenes[agent.species] * model.ploidy[agent.species])
    mutated = true
  end
  # mutate pleiotropy matrix
  if rand(model.rng, model.mutProbs[agent.species][2])
    randnumbers = rand(model.rng, model.rng, model.mutMagnitudes[agent.species][2], size(agent.pleiotropyMat))
    agent.pleiotropyMat[randnumbers] .= .!agent.pleiotropyMat[randnumbers]
    mutated = true
  end
  # mutate epistasis matrix
  if rand(model.rng, model.mutProbs[agent.species][3])
    agent.epistasisMat .+= rand(model.rng, model.mutMagnitudes[agent.species][3], size(agent.epistasisMat))
    mutated = true
  end
  # update biotic and abiotic phenotypes and W
  if mutated
    update_fitness!(agent, model)
  end
end

coord2vertex(c, model) = c[1] + (size(model.space)[1] * (c[2] - 1))

function vertex2coord(vertex::Int, dims::Tuple{Int,Int})
  x = vertex % dims[1]
  if x == 0
    x = dims[1]
  end
  y = ceil(Int, vertex / dims[1])
  return (x, y)
end

vertex2coord(v::Int, model::ABM) = vertex2coord(v, size(model.space))


function update_fitness!(agent::Ind, model::ABM)
  abiotic_phenotype = get_abiotic_phenotype(agent, model)
  biotic_phenotype = get_biotic_phenotype(agent, model)
  W = abiotic_fitness(abiotic_phenotype, agent.species, agent.pos, model)
  W = adjust_fitness(W, agent.species, model)
  agent.biotic_phenotype .= biotic_phenotype
  agent.abiotic_phenotype .= abiotic_phenotype
  agent.W = W
end