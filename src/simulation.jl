
"""
A `struct` for individuals that keeps individual-specific variables.
"""
mutable struct Ind{A<:Integer, B<:AbstractVector, D<:AbstractFloat, E<:AbstractArray} <: AbstractAgent
  id::A  # the individual ID
  species::A  # the species ID the individual belongs to
  y::B  # a vector that specifies the importance of each trait
  W::D  # fitness. W = exp(γ .* transpose(z .- θ)*inv(ω)*(z .- θ)). # z::C  # the phenotype vector. z = By .+ μ
  B::E  # the B matrix
end

"""
    model_initiation(;L, P, B, γ, m, T, Ω, M, MB, N, Y, E, generations, seed=0)

Innitializes the model.

# Parameters and their example values

* L a tuple  specifying the number of loci l for each species
* P a tuple  specifying the number of traits p for each species
* B  [Random.bitrand(i[1], i[2]) for i in zip(P, L)]  a tuple  of pleiotropy matrices, one for  each species. Each matrix consists of zeros and ones only. make sure no rows are all zero (a trait * is not controled by any locus)
* γ [-0.5, -0.5] a tuple  of selection coefficients for each species
* m ploidy, m=2 diploid
* T [randn(Float16, n) for n in P]  a tuple  of arrays, each θ specifying optimal phenotypes for  each species
* Ω [Symmetric(rand(Float16, i[1], i[2])) for i in zip(P, P)] a tuple  of matrices, each of which ω represents a covariance matrix of the selection surface
* M [0.02, 0.02] A tuple of mutation rates μ for each species
* MB [0.05, 0.05] A tuple of mutation rates μ<sub>B</sub> for each species
* N [1000, 1000] a tuple  for population size of each species
* Y [rand(Float16, i*m) for i in L] a tuple  of Arrays, each specifying the initial y vector of  each species
* E [0.8, 0.8] a tuple  of the variance of a normal distribution ε representing environmental * noise for each species.
* generations 100  number of generations to run the simulation
"""
function model_initiation(;L, P, B, γ, m, T, Ω, M, MB, N, Y, E, generations, seed=0)
  if seed >0
    Random.seed!(seed)
  end

  # Some checks for parameters having the correct dimensions
  @assert m < 3  "Ploidy more than 2 is not implemented"
  @assert length.(Y) .% m == Tuple([0 for i in 1:length(Y)]) "Length of elements in Y are not correct. They should a factor of m"
  @assert size.(B, 2) .% m == Tuple([0 for i in 1:length(Y)]) "Length of elements in B are not correct. They should a factor of m"
  @assert length(γ) == length(N) == length(L) == length(B) == length(T) == length(MB) == length(M) == length(Y) == length(E) == length(P) == length(Ω) "L, B, N, γ, T, MB, M, P, Y, Ω, and E should have the same number of elements"

  Ed = [Normal(0.0, i) for i in E]
  # A descrete non-parametric distribtion of uB for each species
  dnps = [DiscreteNonParametric([true, false], [i, 1-i]) for i in MB]
  Mdists = [Normal(0, i) for i in M]

  properties = Dict(:L => L, :P => P, :B => B, :γ => γ, :m => m, :T => T, :Ω => inv.(Ω), :M => Mdists, :MB => dnps, :N => N, :Y => Y, :E => Ed, :generations => generations)
  model = ABM(Ind, properties=properties)
  # create and add agents
  indcounter = 0
  for (sind, n) in enumerate(properties[:N])
    x = properties[:B][sind]*properties[:Y][sind]
    d = properties[:E][sind]
    for ind in 1:n
      indcounter += 1
      z = x .+ rand(d)
      takeabs = abs.(z .- properties[:T][sind])
      # W = exp(γ[sind] * transpose(takeabs)* Ω[sind] *takeabs)
      W = exp(properties[:γ][sind] * transpose(takeabs)*properties[:Ω][sind]*takeabs)
      W = minimum([1e5, W])
      individual = Ind(indcounter, sind, deepcopy(properties[:Y][sind]), W, deepcopy(properties[:B][sind]))
      model.agents[indcounter] = individual
    end
  end

  return model
end

"""
    model_step!(model::ABM)

A function to define what happens within each step of the model.
"""
function model_step!(model::ABM)
  if model.properties[:m] == 2
    sexual_reproduction!(model)
  end
  selection!(model)
end

function agent_step!(agent::Ind, model::ABM)
  mutation!(agent, model)
  update_fitness!(agent, model)
end

function selection!(model::ABM)
  Agents.sample!(model, nagents(model), :W)  # TODO: add a parameter for carrying capacity, so that the population can grow.
end

"""
Choose a random mate and produce 2 offsprings with recombination.
"""
function sexual_reproduction!(model::ABM)
  nparents = nagents(model)
  mates = mate(model)
  for pair in mates
    reproduce!(model.agents[pair[1]], model.agents[pair[2]], model)
  end

  # kill the parents
  for id in 1:nparents
    kill_agent!(model.agents[id], model)
  end
end

"Returns an array of tuples for each pair of agent ids to reproduce"
function mate(model::ABM)
  same_species = [[k for (k, v) in model.agents if v.species == i] for i in 1:length(model.properties[:N])]

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
  child = add_agent!(model, agent1.species, childy, 0.2, childB)  
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