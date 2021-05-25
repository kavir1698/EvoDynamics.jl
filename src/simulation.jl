
"""
A `struct` for individuals that keeps individual-specific variables.
"""
mutable struct Ind{A, B<:AbstractFloat, C<:AbstractArray, D<:AbstractArray, E<:AbstractArray} <: AbstractAgent
  id::Int  # the individual ID
  pos::A  # the individuals position
  species::Int  # the species ID the individual belongs to
  W::B  # fitness. W = exp(γ .* transpose(sum(epistasisMat, dims=2) .- θ)*inv(ω)*(sum(epistasisMat, dims=2) .- θ)).
  epistasisMat::C  # epistasis matrix
  pleiotropyMat::D  # pleiotropy matrix
  q::E  # expression array
end

"""
    model_initiation(;ngenes, nphenotypes, epistasisMat, pleiotropyMat, expressionArrays, selectionCoeffs, ploidy, optPhenotypes, covMat, mutProbs, N, E, growthrates, interactionCoeffs, mutMagnitudes, generations, migration_rates, K, space=nothing, periodic=false, metric=:chebyshev, seed=0)

Innitializes the model.
"""
function model_initiation(;ngenes, nphenotypes, epistasisMat, pleiotropyMat, expressionArrays, selectionCoeffs, ploidy, optPhenotypes, covMat, mutProbs, N, E, growthrates, interactionCoeffs, mutMagnitudes, generations, migration_rates, K, space=nothing, periodic=false, metric=:chebyshev, interaction_equation="lotkaVoltera_competition", seed=0)
  if seed >0
    Random.seed!(seed)
  end

  # TODO: Users should define migration trait and specify its genes with pleiotropy.
  # TODO: We cannot have an optimal phenotype for migration to be selected, right?
  # covMat, however, should include migration trait
  
  if isnothing(space)
    fspace = GridSpace((1, 1))
  elseif typeof(space) <: Tuple
    fspace = GridSpace(space, periodic=periodic, metric=metric)
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
  @assert length(keys(K)) >= nnodes(fspace) "K should have a key for every node"
  for (k, v) in N
    @assert length(v) == nspecies "Each value in N should have size equal to number of species"
  end
  if !isnothing(migration_rates)
    for item in migration_rates
      if typeof(item) <: AbstractArray
        @assert size(item, 1) == nnodes(fspace) "migration_rates has different rows than there are nodes in space."
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
  for i in eachindex(epistasisMat)
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

  epistasisMatS = [MArray{Tuple{size(epistasisMat[i])...}}(newA[i]) for i in eachindex(newA)]
  pleiotropyMatS = [MArray{Tuple{size(pleiotropyMat[i])...}}(pleiotropyMat[i]) for i in eachindex(pleiotropyMat)]
  expressionArraysS = [MArray{Tuple{size(newQ[i])...}}(newQ[i]) for i in eachindex(newQ)]
  properties = Dict(:ngenes => ngenes, :nphenotypes => nphenotypes, :epistasisMat => epistasisMatS, :pleiotropyMat => pleiotropyMatS, :expressionArrays => expressionArraysS, :growthrates => growthrates, :interactionCoeffs => interactionCoeffs, :selectionCoeffs => selectionCoeffs, :ploidy => ploidy, :optPhenotypes => optPhenotypes, :covMat => inv.(newcovMat), :mutProbs => Mdists, :mutMagnitudes => Ddists, :N => N, :E => Ed, :generations => generations, :K => K, :migration_rates => migration_rates, :nspecies => nspecies, :new_N => N, :interaction_equation => interaction_equation)

  postype = typeof(fspace) <: GridSpace ? typeof(size(fspace)) : typeof(1)
  indtype = EvoDynamics.Ind{postype, typeof(0.1), eltype(properties[:epistasisMat]), eltype(properties[:pleiotropyMat]), eltype(properties[:expressionArrays])}
  model = ABM(indtype, fspace, properties=properties)
  
  model.properties[:nodes] = collect(positions(model))

  # create and add agents
  nns = collect(positions(model))
  for (pos, Ns) in properties[:N]
    for (ind2, n) in enumerate(Ns)
      x = properties[:pleiotropyMat][ind2] * (properties[:epistasisMat][ind2] * properties[:expressionArrays][ind2])  # phenotypic values
      d = properties[:E][ind2]
      for ind in 1:n
        z = x .+ rand(d)
        takeabs = abs.(z .- properties[:optPhenotypes][ind2])
        W = exp(-properties[:selectionCoeffs][ind2] * transpose(takeabs)*properties[:covMat][ind2]*takeabs)[1]
        W = minimum([1e5, W])
        add_agent!(nns[pos], model, ind2, W, MArray{Tuple{size(properties[:epistasisMat][ind2])...}}(properties[:epistasisMat][ind2]), MArray{Tuple{size(properties[:pleiotropyMat][ind2])...}}(properties[:pleiotropyMat][ind2]), MVector{length(properties[:expressionArrays][ind2])}(properties[:expressionArrays][ind2]))
      end
    end
  end

  return model
end

nnodes(x::GridSpace) = prod(size(x))
nnodes(x::GraphSpace) = Agents.LightGraphs.nv(x)
nnodes(x::ABM) = nnodes(x.space)

"""
    model_step!(model::ABM)

A function to define what happens within each step of the model.
"""
function model_step!(model::ABM)
  if sum(model.ploidy) > length(model.ploidy) # there is at least one diploid
    sexual_reproduction!(model)
  end
  selection!(model)
  if model.interaction_equation == "lotkaVoltera_generalized"
    for node in 1:nnodes(model)
      model.new_N[node] = Tuple(lotkaVoltera_generalized(model, node))
    end
  end
end

function agent_step!(agent::Ind, model::ABM)
  mutation!(agent, model)
  migration!(agent, model)
end

function selection!(model::ABM)
  for node in 1:nnodes(model)
    for species in 1:model.nspecies
      sample!(model, species, node, :W)
    end
  end
end

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

"""
For sexual reproduction of diploids.

An offspring is created from gametes that include one allele from each loci and the corresponding column of the epistasisMat matrix.
Each gamete is half of `epistasisMat` (column-wise)
"""
function reproduce!(agent1::Ind, agent2::Ind, model::ABM)
  nloci = Int(model.ngenes[agent1.species]/2)
  loci_shuffled = shuffle(1:nloci)
  loci1 = 1:ceil(Int, nloci/2)
  noci1_dip = vcat(loci_shuffled[loci1], loci_shuffled[loci1] .+ nloci)
  childA = MArray{Tuple{size(agent2.epistasisMat)...}}(agent2.epistasisMat)
  childA[:, noci1_dip] .= agent1.epistasisMat[:, noci1_dip]
  childB = MArray{Tuple{size(agent2.pleiotropyMat)...}}(agent2.pleiotropyMat)
  childB[:, noci1_dip] .= agent1.pleiotropyMat[:, noci1_dip]
  childq = MVector{length(agent2.q)}(agent2.q)
  childq[noci1_dip] .= agent1.q[noci1_dip]
  child = add_agent!(agent1.pos, model, agent1.species, 0.2, childA, childB, childq)  
  update_fitness!(child, model)
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
  takeabs = abs.((Fmat .+ rand(d)) .- model.optPhenotypes[agent.species])
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

"Move the agent to a new node with probabilities given in `migration_rates`"
function migration!(agent::Ind, model::ABM)
  if isnothing(model.migration_rates[agent.species])
    return
  end
  vertexpos = coord2vertex(agent.pos, model)  # cell order. Only works for 2D grids
  row = view(model.migration_rates[agent.species], :, vertexpos)
  new_node = sample(1:length(row), Weights(row))
  if new_node != vertexpos
    move_agent!(agent, model.nodes[new_node], model)
  end
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

"""
Calculates the next population size of a species given its current size, intrinsic growth rate, carrying capacity, and competition with other species.

# Arguments

* node: node number
* species: species number

"""
function lotkaVoltera_competition(model::ABM, species::Int, node::Int)
  Ns = nagents_species(model, node)
  N = Ns[species]
  if N == 0
    return
  end
  cc = model.interactionCoeffs
  species_ids = 1:model.nspecies
  if isnothing(cc) || length(species_ids) == 0
    r = model.growthrates[species]
    K = model.K[node][species]
    nextN = N + r*N*(1 - (N/K))
    return round(Int, nextN)
  else 
    ccs = view(cc, :, species)
    aNs = ccs' * Ns
    aNs -= ccs[species] * Ns[species]  # removes the current species from aNs.
    r = model.growthrates[species]
    K = model.K[node][species]
    nextN = N + r*N*(1 - ((N+aNs)/K))
    return round(Int, nextN)
  end
end

"""
Same as lotkaVoltera_competition but uses the generalized equation (see below). This version calculates the next population sizes for all species.

* xt+1 = D(xt) + D(xt)(r + Axt)
* xt: array of n population densities (N/K) at time t.
* r: vector of n intrinsic growth rates, measuring the growth of population i when grown alone
* A: n*n matrix of interaction coefficients (`model.interactionCoeffs`)
* D(x): diagonal matrix with xt on diagonal
"""
function lotkaVoltera_generalized(model::ABM, node::Int)
  x = nagents_species(model, node) ./ model.K[node]
  r = model.growthrates
  Dx = diagm(x)
  A = model.interactionCoeffs
  new_x = x .+ (Dx*(r .+ A*x))
  new_N = new_x .* model.K[node]
  return round.(Int, new_N)
end

"Returns population size per species in the node"
function nagents_species(model::ABM, node::Int)
  counter = zeros(Int, model.nspecies)
  for ag in agents_in_position(node, model)
    counter[ag.species] += 1
  end
  return counter
end