
"""
A `struct` for individuals that keeps individual-specific variables.
"""
mutable struct Ind{B<:AbstractFloat, C<:AbstractArray, D<:AbstractArray, E<:AbstractArray} <: AbstractAgent
  id::Int  # the individual ID
  pos::Tuple{Int, Int}  # the individuals position
  species::Int  # the species ID the individual belongs to
  W::B  # fitness. W = exp(γ .* transpose(sum(epistasisMat, dims=2) .- θ)*inv(ω)*(sum(epistasisMat, dims=2) .- θ)).
  epistasisMat::C  # epistasis matrix
  pleiotropyMat::D  # pleiotropy matrix
  q::E  # expression array
end

"""
    model_initiation(;ngenes, nphenotypes, epistasisMat, pleiotropyMat, expressionArrays, selectionCoeffs, ploidy, optPhenotypes, covMat, mutProbs, N, E, growthrates, competitionCoeffs, mutMagnitudes, generations, migration_rates, K, space=nothing, periodic=false, moore=false, seed=0)

Innitializes the model.
"""
function model_initiation(;ngenes, nphenotypes, epistasisMat, pleiotropyMat, expressionArrays, selectionCoeffs, ploidy, optPhenotypes, covMat, mutProbs, N, E, growthrates, competitionCoeffs, mutMagnitudes, generations, migration_rates, K, space=nothing, periodic=false, moore=false, seed=0)
  if seed >0
    Random.seed!(seed)
  end
  
  if space == nothing
    fspace = GridSpace((1, 1))
  elseif typeof(space) <: NTuple
    fspace = GridSpace(space, periodic=periodic, moore=moore)
  elseif typeof(space) <: AbstractGraph
    fspace = GraphSpace(space)
  end
  nspecies = length(ngenes)

  # Some checks for parameters having the correct dimensions
  for i in ploidy
    @assert i < 3  "Ploidy more than 2 is not implemented"
  end
  for i in size.(epistasisMat, 2) .% ploidy
    @assert i == 0 "number of columns in epistasisMat are not correct. They should a factor of ploidy"
  end
  @assert length(selectionCoeffs) == length(ngenes) == length(epistasisMat) == length(optPhenotypes) == length(mutProbs) == length(E) == length(nphenotypes) == length(covMat) == length(growthrates) == length(mutMagnitudes) "ngenes, epistasisMat, selectionCoeffs, optPhenotypes, mutProbs, nphenotypes, covMat, growthrates, mutMagnitudes and E should have the same number of elements"
  @assert length(keys(K)) >= nv(fspace) "K should have a key for every node"
  for (k, v) in N
    @assert length(v) == nspecies "Each value in N should have size equal to number of species"
  end
  if migration_rates != nothing
    for item in migration_rates
      if typeof(item) <: AbstractArray
        @assert size(item, 1) == nv(fspace) "migration_rates has different rows than there are nodes in space."
      end
    end
  end

  Ed = [Normal(0.0, i) for i in E]
  Mdists = [[DiscreteNonParametric([true, false], [i, 1-i]) for i in arr] for arr in mutProbs]  # μ (probability of change)
  Ddists = [[Normal(0, ar[1]), DiscreteNonParametric([true, false], [ar[2], 1-ar[2]]), Normal(0, ar[3])] for ar in mutMagnitudes]  # amount of change in case of mutation
  
  # make single-element arrays 2D so that linAlg functions will work
  newA = Array{Array{Float64}}(undef, length(epistasisMat))
  newQ = Array{Array{Float64}}(undef, length(epistasisMat))
  newcovMat = Array{Array{Float64}}(undef, length(epistasisMat))
  for i in 1:length(epistasisMat)
    if length(epistasisMat[i]) == 1
      newA[i] = reshape(epistasisMat[i], 1, 1)
      newQ[i] = reshape(expressionArrays[i], 1, 1)
      newcovMat[i] = reshape(covMat[i], 1, 1)
    else
      newA[i] = epistasisMat[i]
      newQ[i] = expressionArrays[i]
      newcovMat[i] = covMat[i]
    end
  end

  properties = Dict(:ngenes => ngenes, :nphenotypes => nphenotypes, :epistasisMat => newA, :pleiotropyMat => pleiotropyMat, :expressionArrays => newQ, :growthrates => growthrates, :competitionCoeffs => competitionCoeffs, :selectionCoeffs => selectionCoeffs, :ploidy => ploidy, :optPhenotypes => optPhenotypes, :covMat => inv.(newcovMat), :mutProbs => Mdists, :mutMagnitudes => Ddists, :N => N, :E => Ed, :generations => generations, :K => K, :migration_rates => migration_rates, :nspecies => nspecies)

  indtype = EvoDynamics.Ind{typeof(0.1), eltype(properties[:epistasisMat]), eltype(properties[:pleiotropyMat]), eltype(properties[:expressionArrays])}
  # EvoDynamics.Ind{Float64,Array{Float64,2},Array{Bool,2},Array{Float64,1}}
  model = ABM(indtype, fspace, properties=properties)
  
  # create and add agents
  for (pos, Ns) in properties[:N]
    for (ind2, n) in enumerate(Ns)
      x = properties[:pleiotropyMat][ind2] * (properties[:epistasisMat][ind2] * properties[:expressionArrays][ind2])  # phenotypic values
      d = properties[:E][ind2]
      for ind in 1:n
        z = x .+ rand(d)
        takeabs = abs.(z .- properties[:optPhenotypes][ind2])
        W = exp(-properties[:selectionCoeffs][ind2] * transpose(takeabs)*properties[:covMat][ind2]*takeabs)[1]
        W = minimum([1e5, W])
        add_agent!(pos, model, ind2, W, MArray{Tuple{size(properties[:epistasisMat][ind2])...}}(properties[:epistasisMat][ind2]), MArray{Tuple{size(properties[:pleiotropyMat][ind2])...}}(properties[:pleiotropyMat][ind2]), MVector{length(properties[:expressionArrays][ind2])}(properties[:expressionArrays][ind2]))
      end
    end
  end

  return model
end

"""
    model_step!(model::ABM)

A function to define what happens within each step of the model.
"""
function model_step!(model::ABM)
  if sum(model.properties[:ploidy]) > length(model.properties[:ploidy]) # there is at least one diploid
    sexual_reproduction!(model)
  end
  selection!(model)
end

function agent_step!(agent::Ind, model::ABM)
  mutation!(agent, model)
  migration!(agent, model)
end

function selection!(model::ABM)
  for node in 1:nv(model)
    for species in 1:model.properties[:nspecies]
      sample!(model, species, node, :W)
    end
  end
end

"""
Choose a random mate and produce one offsprings with recombination.
"""
function sexual_reproduction!(model::ABM, node_number::Int)
  node_content = get_node_contents(node_number, model)
  mates = mate(model, node_number)
  for pair in mates
    reproduce!(model.agents[pair[1]], model.agents[pair[2]], model)
  end

  # kill the parents
  for id in node_content
    kill_agent!(model.agents[id], model)
  end
end

function sexual_reproduction!(model::ABM)
  for node in 1:nv(model)
    sexual_reproduction!(model, node)
  end
end

"Returns an array of tuples for each pair of agent ids to reproduce"
function mate(model::ABM, node_number::Int)
  node_content = get_node_contents(node_number, model)
  same_species = [[k for k in node_content if model.agents[k].species == i] for i in 1:model.properties[:nspecies] if model.properties[:ploidy][i] == 2]

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

"""
For sexual reproduction of diploids.

An offspring is created from gametes that include one allele from each loci and the corresponding column of the epistasisMat matrix.
Each gamete is half of `epistasisMat` (column-wise)
"""
function reproduce!(agent1::Ind, agent2::Ind, model::ABM)
  nloci = Int(model.properties[:ngenes][agent1.species]/2)
  loci_shuffled = shuffle(1:nloci)
  loci1 = 1:ceil(Int, nloci/2)
  noci1_dip = vcat(loci_shuffled[loci1], loci_shuffled[loci1] .+ nloci)
  childA = MVector{length(agent2.epistasisMat)}(agent2.epistasisMat)
  childA[:, noci1_dip] .= agent1.epistasisMat[:, noci1_dip]
  childB = MVector{length(agent2.pleiotropyMat)}(agent2.pleiotropyMat)
  childB[:, noci1_dip] .= agent1.pleiotropyMat[:, noci1_dip]
  childq = MVector{length(agent2.q)}(agent2.q)
  childq[noci1_dip] .= agent1.q[noci1_dip]
  child = add_agent!(agent1.pos, model, agent1.species, 0.2, childA, childB, childq)  
  update_fitness!(child, model)
end

"Mutate an agent."
function mutation!(agent::Ind, model::ABM)
  # mutate gene expression
  if rand(model.properties[:mutProbs][agent.species][1])
    agent.q .+= rand(model.properties[:mutMagnitudes][agent.species][1], model.properties[:ngenes][agent.species])
  end
  # mutate pleiotropy matrix
  if rand(model.properties[:mutProbs][agent.species][2])
    randnumbers = rand(model.properties[:mutMagnitudes][agent.species][2], size(agent.pleiotropyMat))
    agent.pleiotropyMat[randnumbers] .= .!agent.pleiotropyMat[randnumbers]
  end
  # mutate epistasis matrix
  if rand(model.properties[:mutProbs][agent.species][3])
    agent.epistasisMat .+= rand(model.properties[:mutMagnitudes][agent.species][3], size(agent.epistasisMat))
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
  d = model.properties[:E][agent.species]
  Fmat = agent.pleiotropyMat * (agent.epistasisMat * agent.q)
  takeabs = abs.((Fmat .+ rand(d)) .- model.properties[:optPhenotypes][agent.species])
  W = exp(-model.properties[:selectionCoeffs][agent.species] * transpose(takeabs) * model.properties[:covMat][agent.species] * takeabs)[1]
  W = min(1e5, W)
  agent.W = W
end

function update_fitness!(model::ABM)
  for agent in values(model.agents)
    update_fitness!(agent, model)
  end
end

"Move the agent to a new node with probabilities given in :migration_rates"
function migration!(agent::Ind, model::ABM)
  if model.properties[:migration_rates][agent.species] == nothing
    return
  end
  vertexpos = Agents.coord2vertex(agent.pos, model)
  row = model.properties[:migration_rates][agent.species][vertexpos, :]
  new_node = sample(1:length(row), Weights(row))
  if new_node != vertexpos
    move_agent!(agent, new_node, model)
  end
end

"Replace Inds in a node with a weighted sample by replacement of Inds"
function sample!(model::ABM, species::Int, node_number::Int,
  weight=nothing; replace=true,
  rng::AbstractRNG=Random.GLOBAL_RNG)

  max_id = maximum(keys(model.agents))
  node_content_all = get_node_contents(node_number, model)
  node_content = [i for i in node_content_all if model.agents[i].species == species]
  if length(node_content) == 0
    return
  end
  n = lotkaVoltera(model, species, node_number)
  if n == 0
    return
  end
  if weight != nothing
      weights = Weights([getproperty(model.agents[a], weight) for a in node_content])
      newids = sample(rng, node_content, weights, n, replace=replace)
  else
      newids = sample(rng, node_content, n, replace=replace)
  end

  for id in newids # add new agents to the model
    ag = deepcopy(model.agents[id])
    ag.id = max_id + 1
    add_agent_pos!(ag, model)
    # model.agents[max_id + 1] = deepcopy(model.agents[id])
    # push!(model.space.agent_positions[coord2vertex(model.agents[id].pos, model)], max_id + 1)
    max_id += 1
  end

  # kill old/extra agents
  genocide!(model, node_content)
end

"""
  genocide!(model::ABM, n::Array)
Kill the agents of the model whose IDs are in n.
"""
function genocide!(model::ABM, n::AbstractArray)
  for k in n
    kill_agent!(model.agents[k], model)
  end
end

"""
Calculates the next population size of a species given its current size, intrinsic growth rate, carrying capacity, and competition with other species.

# Arguments

* node: node number
* species: species number

"""
function lotkaVoltera(model::ABM, species::Int, node::Int)
  Ns = nagents_species(model, node)
  N = Ns[species]
  if N == 0
    return
  end
  as = model.properties[:competitionCoeffs]
  species_ids = 1:model.properties[:nspecies]
  if as == nothing || length(species_ids) == 0
    r = model.properties[:growthrates][species]
    K = model.properties[:K][node][species]
    nextN = N + r*N*(1 - (N/K))
    return round(Int, nextN)
  else 
    aNs = as[species_ids]' * [Ns[i] for i in species_ids]
    aNs -= as[species] * Ns[species]  # removes the current species from aNs.
    r = model.properties[:growthrates][species]
    K = model.properties[:K][node][species]
    nextN = N + r*N*(1 - ((N+aNs)/K))
    return round(Int, nextN)
  end
end

"Returns population size per species in the node"
function nagents_species(model::ABM, node::Int)
  counter = Dict{Int, Int}()
  for i in 1:model.properties[:nspecies]
    counter[i] = 0
  end
  for id in model.space.agent_positions[node]
    counter[model.agents[id].species] += 1
  end
  return counter
end