
"""
A `struct` for individuals that keeps individual-specific variables.
"""
mutable struct Ind{B<:AbstractFloat, C<:AbstractArray, D<:AbstractArray, E<:AbstractArray} <: AbstractAgent
  id::Int  # the individual ID
  pos::Tuple{Int, Int}  # the individuals position
  species::Int  # the species ID the individual belongs to
  W::B  # fitness. 
  epistasisMat::C  # epistasis matrix
  pleiotropyMat::D  # pleiotropy matrix
  q::E  # expression array
  age::Int
  sex::Bool
  interaction_history::MArray{S, Int} where S # records the last interaction with all species
  energy::B  # determines whether need to feed (energy=0) or not.
end

struct Params{F<:AbstractFloat, I<:Int}
  ngenes::Vector{I}
  nphenotypes::Vector{I}
  epistasisMat::Vector{MArray{S, F, 2, L} where {S<:Tuple, L}}
  pleiotropyMat::Vector{MArray{S, Bool, 2, L} where {S<:Tuple, L}}
  expressionArrays::Vector{MArray{S, Float64, 1, L} where {S<:Tuple, L}}
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
  recombination::Vector{Poisson{F}}
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
  
  properties = create_properties(dd)

  space = dd[:model]["space"]
  if isnothing(space)
    fspace = GridSpace((1, 1))
  elseif typeof(space) <: Tuple
    fspace = GridSpace(space, periodic=dd[:model]["periodic"], metric=dd[:model]["metric"])
  end

  indtype = EvoDynamics.Ind{typeof(0.1), eltype(properties.epistasisMat), eltype(properties.pleiotropyMat), eltype(properties.expressionArrays)}

  model = ABM(indtype, fspace, properties=properties, scheduler=Schedulers.randomly)
  
  # create and add agents
  for sp in 1:model.nspecies
    for (pos, n) in enumerate(model.N[sp])
      for ind in 1:n
        W = abiotic_fitness(sp, pos, model)
        sex = false
        if model.ploidy[sp] == 2
          sex = rand((true, false))
        end
        interaction_history = MVector{model.nspecies, Int}(fill(0, model.nspecies))
        add_agent!(model.nodes[pos], model, sp, W, MArray{Tuple{size(properties.epistasisMat[sp])...}}(properties.epistasisMat[sp]), MArray{Tuple{size(properties.pleiotropyMat[sp])...}}(properties.pleiotropyMat[sp]), MVector{length(properties.expressionArrays[sp])}(properties.expressionArrays[sp]), 0, sex, interaction_history, 0)
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

  properties = Params(ngenes, nphenotypes, epistasisMatS, pleiotropyMatS, expressionArraysS, growthrates, selectionCoeffs, ploidy, optvals, optinds, Mdists, Ddists, Ns, Ed, generations, nspecies, Ns, migration_traits, vision_radius, check_fraction, migration_thresholds, step, nnodes, biotic_phenotyps, abiotic_phenotyps, max_ages, ids, dd[:model]["food sources"], dd[:model]["interactions"], dd[:model]["resources"], recombination)
  
  return properties
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

"""
    model_step!(model::ABM)

A function to define what happens within each step of the model.
# TODO: change the model such that these functions are at agent level.
"""
function model_step!(model::ABM)
  if sum(model.ploidy) > length(model.ploidy) # there is at least one diploid
    sexual_reproduction!(model)
  end
  # selection!(model)
  # if model.interaction_equation == "lotkaVoltera_generalized"
  #   for node in 1:nnodes(model)
  #     model.new_N[node] = Tuple(lotkaVoltera_generalized(model, node))
  #   end
  # end
  
  model.step += 1
end

function agent_step!(agent::Ind, model::ABM)
  mutation!(agent, model)
  migration!(agent, model)
end

# function selection!(model::ABM)
#   for node in 1:nnodes(model)
#     for species in 1:model.nspecies
#       sample!(model, species, node, :W)
#     end
#   end
# end

"""
Choose a random mate and produce one offsprings with recombination.
"""
function sexual_reproduction!(model::ABM, node_number::Int)
  node_content = ids_in_position(node_number, model)
  mates = mate(model, node_content)
  for pair in mates
    reproduce!(model[pair[1]], model[pair[2]], model)
  end

  # kill the parents
  for id in node_content
    kill_agent!(model[id], model)
  end
end

function sexual_reproduction!(model::ABM)
  for node in 1:nnodes(model)
    sexual_reproduction!(model, node)
  end
end

"Returns an array of tuples for each pair of agent ids to reproduce"
function mate(model::ABM, node_content)
  same_species = [[k for k in node_content if model[k].species == i] for i in 1:model.nspecies if model.ploidy[i] == 2]

  mates = Array{Tuple{Int, Int}}(undef, sum(length.(same_species)))

  counter = 1
  for (index, specieslist) in enumerate(same_species)
    for k in specieslist
      m = rand(same_species[index])
      while m == k
        m = rand(same_species[index])
      end
      mates[counter] = (k, m)
      counter += 1
    end
  end
  return mates
end

function same_species(ag1::Ind, ag2::Ind)
  if ag1.species == ag2.species
    return true
  else
    return false
  end
end

"""
Returns a bitarray for sites to be selected from the first (false) and second (true) homologous chromosome.
"""
function crossing_overs(nsites::Int, ncrossing_overs::Int) 
  output = falses(nsites)
  if ncrossing_overs == 0
    return output
  elseif ncrossing_overs ≥ nsites
    output[1:2:end] .= true
    return output
  end
  breaking_points = sample(1:nsites-1, ncrossing_overs, replace=false, ordered=true)
  last = true
  counter = 0
  for site in 1:nsites
    if counter < ncrossing_overs && site == breaking_points[counter+1]
      counter += 1
      output[site] = last
      last = !last 
    else
      output[site] = last
    end
  end
  return output
end

"""
Returns gamets for epistasisMat, pleiotropyMat, and q.

A ametes includes `cross_overs` sites from one homologous chr and the rest from another.the corresponding column of the `epistasisMat` and `pleiotropyMat` matrices.
Each gamete is half of `epistasisMat` and `pleiotropyMat` (column-wise).
"""
function create_gamete(agent, cross_overs, nsites, first::Bool)
  indices1 = 1:nsites
  indices2 = indices1 .+ 7
  if !first
    indices1, indices2 = indices2, indices1
  end
  epistasisMat_gamet = agent.epistasisMat[:, indices1]
  epistasisMat_gamet[:, cross_overs] .= agent.epistasisMat[:, indices2][:, cross_overs]

  pleiotropyMat_gamet = agent.pleiotropyMat[:, indices1]
  pleiotropyMat_gamet[:, cross_overs] .= agent.pleiotropyMat[:, indices2][:, cross_overs]

  # q_gamet = MVector{nsites}(agent.q[indices1])  MArray{Tuple{nsites*2, nsites}}(
  q_gamet = agent.q[indices1]
  q_gamet[cross_overs] = agent.q[indices2][cross_overs]
  
  return epistasisMat_gamet, pleiotropyMat_gamet, q_gamet
end

"""
Adds new individual(s) to the model as offsprings of `ag1` and `ag2`.
# TODO: calculate the number of offsprings and do this for each offspring
"""
function sexual_reproduction(ag1::Int, ag2::Ind, model::ABM)
  species = ag1.species
  nsites = model.ngenes[species]
  if model.recombination[species] == 0
    nco1, nco2 = 0, 0
  else
    nco1, nco2 = rand(model.recombination[species], 2)
  end
  cross_overs1 = crossing_overs(nsites, nco1)
  cross_overs2 = crossing_overs(nsites, nco2)
  gamets1 = create_gamete(ag1, cross_overs1, nsites, rand((true, false)))
  gamets2 = create_gamete(ag2, cross_overs2, nsites, rand((true, false)))

  sex = rand((true, false))
  interaction_history = MVector{model.nspecies, Int}(fill(0, model.nspecies))
  W = abiotic_fitness(species, ag1.pos, model)
  add_agent!(ag1.pos, model, ag1.species, W, hcat(gamets1[1], gamets2[1]), hcat(gamets1[2], gamets2[2]), vcat(gamets1[3], gamets2[3]), 0, sex, interaction_history, 0)
end

function abiotic_fitness(species, pos, model)
  x = model.pleiotropyMat[species] * (model.epistasisMat[species] * model.expressionArrays[species])  # phenotypic values
  d = model.E[species]
  z = x .+ rand(d)
  abp = model.abiotic_phenotypes[species]
  W = abiotic_distance(z[abp], return_opt_phenotype(species, model.step[1], pos, model), variance)
  return W
end

"Mutate an agent."
function mutation!(agent::Ind, model::ABM)
  # mutate gene expression
  if rand(model.mutProbs[agent.species][1])
    agent.q .+= rand(model.mutMagnitudes[agent.species][1], model.ngenes[agent.species])
  end
  # mutate pleiotropy matrix
  if rand(model.mutProbs[agent.species][2])
    randnumbers = rand(model.mutMagnitudes[agent.species][2], size(agent.pleiotropyMat))
    agent.pleiotropyMat[randnumbers] .= .!agent.pleiotropyMat[randnumbers]
  end
  # mutate epistasis matrix
  if rand(model.mutProbs[agent.species][3])
    agent.epistasisMat .+= rand(model.mutMagnitudes[agent.species][3], size(agent.epistasisMat))
  end
  update_fitness!(agent, model)
end

function mutation!(model::ABM)
  for agent in values(model.agents)
    mutation!(agent, model)
  end
end

"Recalculate the fitness of `agent`"
function update_fitness!(agent::Ind, model::ABM)
  d = model.E[agent.species]
  Fmat = agent.pleiotropyMat * (agent.epistasisMat * agent.q)
  takeabs = abs.((Fmat .+ rand(d)) .- return_opt_phenotype(agent.species, model.step, model.optPhenotypes))
  W = exp(-model.selectionCoeffs[agent.species] * transpose(takeabs) * model.covMat[agent.species] * takeabs)[1]
  W = min(1e5, W)
  agent.W = W
end

function update_fitness!(model::ABM)
  for agent in values(model.agents)
    update_fitness!(agent, model)
  end
end

coord2vertex(c, model) = c[1] + (size(model.space)[1] * (c[2]-1))

function migration!(agent::Ind, model::ABM)
  if isnothing(model.migration_traits[agent.species])
    return
  elseif model.migration_thresholds[agent.species] > get_migration_trait(agent, model)
    return
  end

  sites = collect(nearby_positions(agent, model, model.vision_radius[agent.species]))
  nsites = length(sites)
  nsites_selected = round(Int, model.check_fraction[agent.species] * nsites)
  sites_checked = EvoDynamics.sample(sites, nsites_selected, replace=false)
  # TODO: check_site
  # vertexpos = coord2vertex(agent.pos, model)  # cell order. Only works for 2D grids
  # row = view(model.physical_distance[agent.species], :, vertexpos)
  # new_node = sample(1:length(row), Weights(row))
  # TODO: add migration cost. If it survives migration, then do the migration.
  move_agent!(agent, model.nodes[new_node], model)
end

# TODO
"Agent evaluates the site and gives it a score"
function check_site(agent, site, model)
end

function get_migration_trait(agent, model)
  mtrait = agent.pleiotropyMat[model.migration_traits[agent.species], :]
  sum(agent.epistasisMat[mtrait, :] * agent.q)
end

"Replace Inds in a node with a weighted sample by replacement of Inds"
function sample!(model::ABM, species::Int, node_number::Int,
  weight=nothing; replace=true,
  rng::AbstractRNG=Random.GLOBAL_RNG)

  node_content_all = ids_in_position(node_number, model)
  node_content = [i for i in node_content_all if model[i].species == species]
  length(node_content) == 0 && return

  if model.interaction_equation == "lotkaVoltera_competition"
    n = lotkaVoltera_competition(model, species, node_number)
  elseif model.interaction_equation == "lotkaVoltera_generalized"
    n = model.new_N[node_number][species]
  end
  n <= 0 && return  # n can be <0 in lotkaVoltera_generalized

  if !isnothing(weight)
      weights = Weights([getproperty(model[a], weight) for a in node_content])
      newids = sample(rng, node_content, weights, n, replace=replace)
  else
      newids = sample(rng, node_content, n, replace=replace)
  end

  add_newids!(model, node_content, newids)
end

"Used in sample!"
function add_newids!(model, node_content, newids)
  n = nextid(model)
  for id in node_content
    if !in(id, newids)
      kill_agent!(model[id], model)
    else
      noccurances = count(x->x==id, newids)
      for t in 2:noccurances
        newagent = deepcopy(model[id])
        newagent.id = n
        add_agent_pos!(newagent, model)
        n += 1
      end
    end
  end
end

"""
  genocide!(model::ABM, n::Array)
Kill the agents of the model whose IDs are in n.
"""
function genocide!(model::ABM, n::AbstractArray)
  for k in n
    kill_agent!(model[k], model)
  end
end
