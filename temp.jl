
"""
A `struct` for individuals that keeps individual-specific variables.
"""
mutable struct Ind{A<:Integer, B<:AbstractVector, D<:AbstractFloat, E<:AbstractArray} <: AbstractAgent
  id::A  # the individual ID
  pos  # the individuals position
  species::A  # the species ID the individual belongs to
  y::B  # a vector that specifies the importance of each trait
  W::D  # fitness. W = exp(γ .* transpose(z .- θ)*inv(ω)*(z .- θ)). # z::C  # the phenotype vector. z = By .+ μ
  B::E  # the B matrix
end

"""
    model_initiation(;L, P, B, γ, m, T, Ω, M, MB, N, Y, E, generations, seed=0)

Innitializes the model.
"""
function model_initiation(;L, P, B, γ, m, T, Ω, M, MB, N, Y, E, R, C, generations, migration_rates, K, space=nothing, periodic=false, moore=false, seed=0)
  if seed >0
    Random.seed!(seed)
  end
  
  if space == nothing
    fspace = Space((1, 1))
  elseif typeof(space) <: NTuple
    fspace = Space(space, periodic=periodic, moore=moore)
  elseif typeof(space) <: AbstractGraph
    fspace = Space(space)
  end
  nspecies = length(L)

  # Some checks for parameters having the correct dimensions
  for i in m
    @assert i < 3  "Ploidy more than 2 is not implemented"
  end
  @assert length.(Y) .% m == Tuple([0 for i in 1:length(Y)]) "Length of elements in Y are not correct. They should a factor of m"
  @assert size.(B, 2) .% m == Tuple([0 for i in 1:length(Y)]) "Length of elements in B are not correct. They should a factor of m"
  @assert length(γ) == length(L) == length(B) == length(T) == length(MB) == length(M) == length(Y) == length(E) == length(P) == length(Ω) == length(R) "L, B, γ, T, MB, M, P, Y, Ω, R and E should have the same number of elements"
  @assert length(keys(K)) >= nv(fspace) "K should have a key for every node"
  for (k, v) in N
    @assert length(v) == nspecies "Each value in N should have size equal to number of species"
  end
  for item in migration_rates
    if typeof(item) <: AbstractArray
      @assert size(item, 1) == nv(fspace) "migration_rates has different rows than there are nodes in space."
    end
  end

  Ed = [Normal(0.0, i) for i in E]
  # A descrete non-parametric distribtion of uB for each species
  dnps = [DiscreteNonParametric([true, false], [i, 1-i]) for i in MB]
  Mdists = [Normal(0, δ) for δ in M]  # TODO: δ (amount of change) should be different from μ (probability of change).
  
  properties = Dict(:L => L, :P => P, :B => B, :R => R, :C => C, :γ => γ, :m => m, :T => T, :Ω => inv.(Ω), :M => Mdists, :MB => dnps, :N => N, :Y => Y, :E => Ed, :generations => generations, :K => K, :migration_rates => migration_rates, :nspecies => nspecies)
  model = ABM(Ind, fspace, properties=properties)
  # create and add agents
  for (pos, Ns) in properties[:N]
    for (sind, n) in enumerate(Ns)
      x = properties[:B][sind]*properties[:Y][sind]
      d = properties[:E][sind]
      for ind in 1:n
        z = x .+ rand(d)
        takeabs = abs.(z .- properties[:T][sind])
        W = exp(properties[:γ][sind] * transpose(takeabs)*properties[:Ω][sind]*takeabs)
        W = minimum([1e5, W])
        add_agent!(pos, model, sind, deepcopy(properties[:Y][sind]), W, deepcopy(properties[:B][sind]))
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
  if sum(model.properties[:m]) > length(model.properties[:m]) # there is at least one diploid
    sexual_reproduction!(model)
  end
  selection!(model)
end

function agent_step!(agent::Ind, model::ABM)
  mutation!(agent, model)
  update_fitness!(agent, model)
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
  same_species = [[k for k in node_content if model.agents[k].species == i] for i in 1:model.properties[:nspecies] if model.properties[:m][i] == 2]

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

An offspring is created from gametes that include one allele from each loci and the corresponding column of the B matrix.
Each gamete is half of `y` and half of `B` (column-wise)
"""
function reproduce!(agent1::Ind, agent2::Ind, model::ABM)
  nloci = Int(model.properties[:L][agent1.species]/2)
  loci_shuffled = shuffle(1:nloci)
  loci1 = 1:ceil(Int, nloci/2)
  noci1_dip = vcat(loci_shuffled[loci1], loci_shuffled[loci1] .+ nloci)
  childB = deepcopy(agent2.B)
  childy = deepcopy(agent2.y)
  childB[:, noci1_dip] .= agent1.B[:, noci1_dip]
  childy[noci1_dip] .= agent1.y[noci1_dip]
  child = add_agent!(agent1.pos, model, agent1.species, childy, 0.2, childB)  
  update_fitness!(child, model)
end

"Mutate an agent."
function mutation!(agent::Ind, model::ABM)
  # u
  agent.y .+= rand(model.properties[:M][agent.species], model.properties[:L][agent.species]);
  # uB
  randnumbers = rand(model.properties[:MB][agent.species], size(agent.B))
  agent.B[randnumbers] .= .!agent.B[randnumbers]
end

function mutation!(model::ABM)
  for agent in values(model.agents)
    mutation!(agent, model)
  end
end

"Recalculate the fitness of `agent`"
function update_fitness!(agent::Ind, model::ABM)
  d = model.properties[:E][agent.species]
  takeabs = abs.((agent.B * agent.y .+ rand(d)) .- model.properties[:T][agent.species])
  W = exp(model.properties[:γ][agent.species] * transpose(takeabs) * model.properties[:Ω][agent.species] * takeabs)
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
  vertexpos = coord2vertex(agent.pos, model)
  row = model.properties[:migration_rates][agent.species][vertexpos, :]
  new_node = sample(1:length(row), Weights(row))
  if new_node != vertexpos
    move_agent!(agent, new_node, model)
  end
end

"Replace Inds in a node with a weighted by :W sample of Inds with replacement"
function sample!(model::ABM, species::Int, node_number::Int,
  weight=nothing; replace=true,
  rng::AbstractRNG=Random.GLOBAL_RNG)

  max_id = maximum(keys(model.agents))
  node_content_all = get_node_contents(node_number, model)
  node_content = [i for i in node_content_all if model.agents[i].species == species]
  n = lotkaVoltera(model, species, node_number)
  if n == 0

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
  as = model.properties[:C]
  species_ids = collect(1:model.properties[:nspecies])
  splice!(species_ids, species)
  aNs = as[species_ids]' * [Ns[i] for i in species_ids]
  r = model.properties[:R][species]
  K = model.properties[:K][node][species]
  nextN = N + r*N*(1 - ((N+aNs)/K))
  return round(Int, nextN)
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