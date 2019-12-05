# Define model properties
P = (4, 5)
L = (7, 8)
m = 1
model_properties = Dict(
  :L => L, # a tuple  specifying the number of loci l for each species
  :P => P, # a tuple  specifying the number of traits p for each species
  :B =>  Tuple([Random.bitrand(i[1], i[2]) for i in zip(P, L)]),  # a tuple  of pleiotropy matrices, one for each species. Each matrix consists of zeros and ones only. make sure no rows are all zero (a trait is not controled :by any locus)
  :γ => (-0.5, -0.5), # a tuple  of selection coefficients for each species
  :m => m, # ploidy, m=2 diploid. Currently only haploids are implemented
  :T => Tuple([randn(Float16, n) for n in P]), # a tuple  of arrays, each θ specifying optimal phenotypes for each species
  :Ω => Tuple([Symmetric(rand(Float16, i[1], i[2])) for i in zip(P, P)]), # a tuple  of matrices, each of wich ω represents a covariance matrix of the selection surface
  :M => (0.02, 0.02), # a tuple of mutation rates μ for each species
  :MB => (0.05, 0.05), # a tuple of mutation rates μ<sub>B</sub> for each species
  :N => (1000, 1000), # a tuple  for population size of each species
  :Y => Tuple([rand(Float16, i*m) for i in L]), # a tuple  of Arrays, each specifying the initial y vector of each species
  :E => (0.8, 0.8), # a tuple  of the variance of a normal distribution ε representing environmental noise for each species.
  :generations => 100 # number of generations to run the simulation
)

# Data collection params
function mean_fitness(model::ABM)
  nspecies = length(model.properties[:P])
  mean_fitness = Array{Float32}(undef, nspecies)
  for species in 1:nspecies
    fitness = mean([i.W for i in values(model.agents) if i.species == species])
    mean_fitness[species] = fitness
  end

  return Tuple(mean_fitness)
end

collect_properties = Dict(:model => [mean_fitness])
when = 1:model_properties[:generations];

# create model
model = model_initiation(;model_properties...)

# run model and collect data
data = step!(model, dummystep, model_step!, model_properties[:generations], collect_properties, when=when)

for species in 1:length(model.properties[:P])
  data[!, Symbol("fitness_sp_$species")] = getfield.(data[!, 1], species)
end

Agents.select!(data, Agents.Not(Symbol("mean_fitness(model)")))