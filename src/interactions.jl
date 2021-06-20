
"""
Calculates the average distance between two sets of phenotypes.

Distance is a probability and is calculated by the number of standard deviation that phenotype 2 is different from phenotype 1.

Variance can stand for selection coefficient, with a larger variance, less sever the difference.
"""
function phenotypic_distance(phenotype1, phenotype2, sign::Symbol, variance::Float64)
  distance = 0.0
  counter = 0.0
  for ph1 in phenotype1
    for ph2 in phenotype2
      d= abs(0.5 - cdf(Normal(ph1,  variance), ph2)) # 0.5 is the max possible value
      distance += d
      counter += 1.0
    end
  end
  distance += distance  # Double the distance so that it is in range 0 to 1.
  distance /= counter  # average distance
  if sign != :match
    distance = 1.0 - distance
  end
  return distance
end

#= testing the distance function
traits1 = [0.5, 1.5]
traits2 = [0.7, 2.0, 1.3]
variance = 1.0

d = phenotypic_distance(traits1, traits2, :match, variance)
dm = phenotypic_distance(traits1, traits2, :mismatch, variance)
=#

"""
Abiotic fitness
"""
function abiotic_distance(phenotype, optimal, variance)
  distance = 0.0
  counter = 0.0
  for i in eachindex(phenotype)
    d = abs(0.5 - cdf(Normal(optimal[i],  variance), phenotype[i]))
    distance += d
    counter += 1.0
  end
  distance += distance
  distance /= counter
  return distance
end

#= testing
ind1 = [1.5, 1.2]
opt1 = [0.5, 1.0]
variance = 1
abiotic_distance(ind1, opt1, variance)
=#


"""
Change an individuals survivability or decide its number of children.

Each individual has a "raw" survival prob. that is determined by its abiot phenotypes.
Interactions increase or decrease that probablity.
"""
function interaction!(individual1, individual2, model)
  interaction_prob = phenotypic_distance(individual1.biotic, individual2.biotic, model.sign, model.variance)
  ind1abiotic = abiotic_distance(individual1.abiotic, individual1.opt, model.variance)
  ind2abiotic = abiotic_distance(individual2.abiotic, individual1.opt, model.variance)
  if model.inx_type == :competition
    # decrease survival probability of both individuals
    survival_prob1 = ind1abiotic - interaction_prob
    survival_prob2 = ind2abiotic - interaction_prob
  elseif model.inx_type == :cooperation || model.inx_type == :reproduction
    # increase survival probability  of both individuals
    survival_prob1 = ind1abiotic + interaction_prob
    survival_prob2 = ind2abiotic + interaction_prob
  elseif model.inx_type == :predation
    # if successful the predator (individual1) survives and the prey (individual2) dies
    if rand() < interaction_prob
      # individual1 recharges
      # kill!(individual2)
    end
  elseif model.inx_type == :prey
    # inverse of the previous status
    if rand() < interaction_prob
      # individual2 recharges
      # kill!(individual1)
    end
  end

  survival_prob1 *= model.selection_coefficient
  survival_prob2 *= model.selection_coefficient

  return survival_prob1, survival_prob2
end

#= testing

individual1 = (biotic = [0.5, 1.2], abiotic = [2.4, 0.3], opt = [1.5, 0.6])
individual2 = (biotic = [0.7, 2.0, 1.3], abiotic = [0.4, 1.3], opt = [0.5, 1.0])
model = (sign = :match, variance = 1.0, inx_type = :competition, selection_coefficient=0.5, growth_rate=1.2)

s1, s2 = interaction!(individual1, individual2, model)

# if reproduction
n_offsprings = rand(Poisson(mean([s1, s2]) * model.growth_rate))
=#