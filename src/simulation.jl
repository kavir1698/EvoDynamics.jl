
"""
A `struct` for individuals that keeps individual-specific variables.
"""
mutable struct Ind{A<:Integer, B<:AbstractVector, C<:AbstractVector, D<:AbstractFloat, E<:AbstractArray} <: AbstractAgent
  id::A  # the individual ID
  species::A  # the species ID the individual belongs to
  y::B  # a vector that specifies the importance of each trait
  z::C  # the phenotype vector. z = By .+ μ
  W::D  # fitness. W = exp(γ .* transpose(z .- θ)*inv(ω)*(z .- θ))
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
* m ploidy, m=2 diploid. Currently only haploids are implemented
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

  Ed = [Normal(0.0, i) for i in E]

  # A descrete non-parametric distribtion of uB for each species
  dnps = [DiscreteNonParametric([true, false], [i, 1-i]) for i in MB]

  properties = Dict(:L => L, :P => P, :B => B, :γ => γ, :m => m, :T => T, :Ω => inv.(Ω), :M => M, :MB => dnps, :N => N, :Y => Y, :E => Ed, :generations => generations)
  model = ABM(Ind, properties=properties, scheduler=random_activation)
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
      individual = Ind(indcounter, sind, deepcopy(properties[:Y][sind]), z, W, deepcopy(properties[:B][sind]))
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
  selection!(model)
  dists = [Normal(0, i) for i in model.properties[:M]]
  mutation!(model, dists)
end

function selection!(model::ABM)
  Agents.sample!(model, nagents(model), :W)
end

"Mutate all agents."
function mutation!(model::ABM, dists)
  for agent in values(model.agents)
    mutation!(agent, model, dists)
  end
  # recalculate fitness
  update_fitness!(model)
end

"Mutate an agent."
function mutation!(agent::Ind, model::ABM, dists)
  # u
  agent.y .+= rand(dists[agent.species], model.properties[:L][agent.species]);
  # uB
  randnumbers = rand(model.properties[:MB][agent.species], size(agent.B))
  agent.B[randnumbers] .= .!agent.B[randnumbers]
end

function update_fitness!(model::ABM)
  for agent in values(model.agents)
    update_fitness!(agent, model)
  end
end

function update_fitness!(agent::Ind, model::ABM)
  d = model.properties[:E][agent.species]
  takeabs = (agent.B * agent.y .+ rand(d)) .- model.properties[:T][agent.species]
  W = exp(model.properties[:γ][agent.species] * transpose(takeabs) * model.properties[:Ω][agent.species] * takeabs)
  W = min(1e5, W)
  agent.W = W
end
