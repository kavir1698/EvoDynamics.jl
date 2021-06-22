
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
      if ag1.sex != ag2.sex && model.ploidy[sp1] == 2 # reproduce
        reproduce!(ag1::Ind, ag2::Ind, model::ABM)
      else # change in fitness
        inx_prob = interaction_power(ag1, ag2, model)
        abiotic1 = abiotic_fitness(ag1, model)
        abiotic2 = abiotic_fitness(ag2, model)
        ix_dir = sign(ix_value)
        
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
  breaking_points = sample(1:(nsites-1), ncrossing_overs, replace=false, ordered=true)
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
Returns gametes for epistasisMat, pleiotropyMat, and q.

A ametes includes `cross_overs` sites from one homologous chr and the rest from another.the corresponding column of the `epistasisMat` and `pleiotropyMat` matrices.
Each gamete is half of `epistasisMat` and `pleiotropyMat` (column-wise).
"""
function create_gamete(agent, cross_overs, nsites, first::Bool)
  indices1 = 1:nsites
  indices2 = indices1 .+ nsites
  if !first
    indices1, indices2 = indices2, indices1
  end
  epistasisMat_gamete = agent.epistasisMat[:, indices1]
  epistasisMat_gamete[:, cross_overs] .= agent.epistasisMat[:, indices2][:, cross_overs]

  pleiotropyMat_gamete = agent.pleiotropyMat[:, indices1]
  pleiotropyMat_gamete[:, cross_overs] .= agent.pleiotropyMat[:, indices2][:, cross_overs]

  q_gamete = agent.q[indices1]
  q_gamete[cross_overs] .= agent.q[indices2][cross_overs]
  
  return epistasisMat_gamete, pleiotropyMat_gamete, q_gamete
end

"""
Adds new individual(s) to the model as offsprings of `ag1` and `ag2`.
"""
function create_one_offspring(ag1::Ind, ag2::Ind, model::ABM)
  species = ag1.species
  nsites = model.ngenes[species]

  if model.recombination[species] == 0
    nco1, nco2 = 0, 0
  else
    nco1, nco2 = rand(model.recombination[species], 2)
  end

  cross_overs1 = crossing_overs(nsites, nco1)
  cross_overs2 = crossing_overs(nsites, nco2)
  gametes1 = create_gamete(ag1, cross_overs1, nsites, rand((true, false)))
  gametes2 = create_gamete(ag2, cross_overs2, nsites, rand((true, false)))

  sex = rand((true, false))
  interaction_history = MVector{model.nspecies, Int}(fill(0, model.nspecies))
  
  # Merge gametes
  nsites2 = nsites*2
  episMat = MMatrix{nsites2, nsites2}(hcat(gametes1[1], gametes2[1]))
  pleioMat = MMatrix{model.nphenotypes[species],nsites2}(hcat(gametes1[2], gametes2[2]))
  q = MVector{nsites2}(vcat(gametes1[3], gametes2[3]))
  abph = get_abiotic_phenotype(species, episMat, pleioMat, q, model) 
  bph = get_biotic_phenotype(species, episMat, pleioMat, q, model)
  W = abiotic_fitness(abph, species, ag1.pos, model)
  initial_energy = model.initial_energy[species]

  add_agent!(ag1.pos, model, ag1.species, bph, abph, episMat, pleioMat, q, 0, sex, interaction_history, initial_energy, W)
end

function reproduce!(ag1::Ind, ag2::Ind, model::ABM)
  reproduction_success = 1.0 - phenotypic_distance(ag1, ag2, model)
  growth_rate = model.growthrates[ag1.species]
  nchildren = rand(Poisson(reproduction_success * growth_rate))
  for c in 1:nchildren
    create_one_offspring(ag1, ag2, model)
  end
end
