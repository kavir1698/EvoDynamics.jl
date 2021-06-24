
function get_biotic_phenotype(species, epistasisMat, pleiotropyMat, expressionArray, model::ABM)
  abp = model.biotic_phenotypes[species]
  x = pleiotropyMat[abp, abp] * (epistasisMat[abp, abp] * expressionArray[abp])  # phenotypic values
  d = model.E[species]
  z = x .+ rand(d)
  return z
end

get_biotic_phenotype(ag::Ind, model::ABM) = get_biotic_phenotype(ag.species, ag.epistasisMat, ag.pleiotropyMat, ag.q, model)

function get_abiotic_phenotype(species, epistasisMat, pleiotropyMat, expressionArray, model::ABM)
  abp = model.abiotic_phenotypes[species]
  x = pleiotropyMat[abp, abp] * (epistasisMat[abp, abp] * expressionArray[abp])  # phenotypic values
  d = model.E[species]
  z = x .+ rand(d)
  return z
end

get_abiotic_phenotype(ag::Ind, model::ABM) = get_abiotic_phenotype(ag.species, ag.epistasisMat, ag.pleiotropyMat, ag.q, model)


"""
Calculates the average distance between two sets of phenotypes.

Distance is a probability and is calculated by the number of standard deviation that phenotype 2 is different from phenotype 1.

Variance can stand for selection coefficient, with a larger variance, less sever the difference.
"""
function phenotypic_distance(phenotype1, phenotype2, interaction_coeff, variance::Real)
  distance = 0.0
  counter = 0.0
  for ph1 in phenotype1
    for ph2 in phenotype2
      d= abs(0.5 - cdf(Normal(ph1,  variance), ph2)) # 0.5 is the value when ph1 and ph2 are identical
      distance += d
      counter += 1.0
    end
  end
  distance += distance  # Double the distance so that it is in range 0 to 1.
  distance /= counter  # average distance
  if interaction_coeff < 0
    distance = 1.0 - distance
  end
  return distance
end

phenotypic_distance(ag1::Ind, ag2::Ind, model) = phenotypic_distance(ag1.biotic_phenotype, ag2.biotic_phenotype, model.interactions[ag1.species, ag2.species], variance)

"""
Abiotic distance to the optimal values
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

function abiotic_fitness(abiotic_phenotype, species::Int, pos, model::ABM)
  rawW = 1.0 - abiotic_distance(abiotic_phenotype, return_opt_phenotype(species, model.step[1], pos, model), variance)
  W = adjust_abiotic_fitness(rawW, model.selectionCoeffs[species])
  return W
end

abiotic_fitness(ag::Ind, model::ABM) = abiotic_fitness(ag.abiotic_phenotype, ag.species, ag.pos, model)

adjust_abiotic_fitness(rawfitness::AbstractFloat, selection_coeff::AbstractFloat) = 1.0 - ((1.0 - rawfitness) * selection_coeff)

function interaction_power(ag1, ag2, model)
  d = phenotypic_distance(ag1, ag2, model)
  prob = (1.0 - d) * abs(model.interactions[ag1.species, ag2.species])
  return prob
end

"""
    interact!(ag1::Ind, ag2::Ind, model::ABM)

Changes the two interacting individuals' survivability or decide their number of children.

Each individual has a "raw" survival prob. that is determined by its abiot phenotypes.
Interactions increase or decrease that probablity.
"""
function interact!(ag1::Ind, ag2::Ind, model::ABM)
  sp1 = ag1.species
  sp2 = ag2.species

  # record interaction history
  ag1.interaction_history[sp2] = model.step[1]
  ag2.interaction_history[sp1] = model.step[1]

  if sp1 != sp2 && model.food_sources[sp1, sp2] > 0  # predation
    d = phenotypic_distance(ag1, ag2, model)
    if rand() < d
      eat!(ag1, ag2, model)
    end
  elseif sp1 != sp2 && model.food_sources[sp2, sp1] > 0 # predation
    d = phenotypic_distance(ag1, ag2, model)
    if rand() < d
      eat!(ag2, ag1, model)
    end
  else # interaction
    ix_value = model.interactions[sp1, sp2]
    if ix_value != 0
      if sp1 == sp2 && ag1.sex != ag2.sex && model.ploidy[sp1] == 2 # reproduce
        reproduce!(ag1::Ind, ag2::Ind, model::ABM)
      else # change in fitness
        inx_prob = interaction_power(ag1, ag2, model)
        abiotic1 = abiotic_fitness(ag1, model)
        abiotic2 = abiotic_fitness(ag2, model)
        ix_dir = sign(ix_value)
        
        # update agents' fitness
        ag1.W = abiotic1 + (inx_prob * ix_dir)
        ag2.W = abiotic2 + (inx_prob * ix_dir)
      end
    end
  end
end

"""
`ag1` eats `ag2`.
"""
function eat!(ag1, ag2, model)
  ag1.energy += 1
  kill_agent!(ag2, model)
end

"""
Returns as many random ids from the agent site as the number of species.
"""
function target_species_ids(agent, model::ABM)
  nspecies = model.nspecies
  allids = collect(nearby_ids(agent, model, 0))
  foundlen = length(allids)
  if foundlen == 0
    return Int[]
  elseif foundlen < nspecies
    return allids
  else
    species_ids = sample(allids, nspecies, replace=false)
    return species_ids  
  end
end

"""
    interact!(agent::Ind, model::ABM)

Includes interaction with other species such as competition, cooperation, and predation. It also includes sexual reproduction (if meeting an individual from other sex) and within species competition/cooperation (if meeting an individual from the same sex).
"""
function interact!(agent::Ind, model::ABM)
  sp = agent.species
  target_ids = target_species_ids(agent, model)
  for id in target_ids
    target = model[id]
    target_sp = target.species
    if agent.interaction_history[target_sp] != model.step[1] || target.interaction_history[sp] != model.step[1] # if agent and target have not interacted with such species before
      interact!(agent, target, model)
      # if agent was a prey, check whether it is still alive
      if model.food_sources[target_sp, sp] > 0
        if !haskey(model.agents, agent.id)
          return
        end
      end
    end
  end
end