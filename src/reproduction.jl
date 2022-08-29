
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

A gametes includes `cross_overs` sites from one homologous chr and the rest from another corresponding column of the `epistasisMat` and `pleiotropyMat` matrices.
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
  interaction_history = MVector{model.nspecies,Int}(fill(0, model.nspecies))

  # Merge gametes
  nsites2 = nsites * 2
  episMat = MMatrix{nsites2,nsites2}(hcat(gametes1[1], gametes2[1]))
  pleioMat = MMatrix{model.nphenotypes[species],nsites2}(hcat(gametes1[2], gametes2[2]))
  q = MVector{nsites2}(vcat(gametes1[3], gametes2[3]))
  abph = get_abiotic_phenotype(species, episMat, pleioMat, q, model)
  bph = get_biotic_phenotype(species, episMat, pleioMat, q, model)
  W = abiotic_fitness(abph, species, ag1.pos, model)
  W = adjust_fitness(W, ag1.species, model)
  initial_energy = model.initial_energy[species]

  offspring = add_agent!(ag1.pos, model, ag1.species, bph, abph, episMat, pleioMat, q, 0, sex, interaction_history, initial_energy, W, true)
  return offspring
end

"""
Sexual reproduction for diploid individuals.
Whether it is assortative or dissortative depends on the sign of the corresponding element in the interaction matrix. 
"""
function reproduce!(ag1::Ind, ag2::Ind, model::ABM)
  if model.mating_schemes[ag1.species] == 0.0
    reproduction_success = 1.0
  else
    d = phenotypic_distance(ag1.biotic_phenotype, ag2.biotic_phenotype, model.mating_schemes[ag1.species], model.biotic_variances[ag1.species])
    reproduction_success = interaction_power(ag1, ag2, d, model)
  end
  growth_rate = model.growthrates[ag1.species]
  nchildren = rand(Poisson(reproduction_success * growth_rate))
  for c in 1:nchildren
    offspring = create_one_offspring(ag1, ag2, model)
    mutate!(offspring, model)
  end
end

"""
Asexual reproduction for the haploid
"""
function reproduce!(agent::Ind, model::ABM)
  if model.ploidy[agent.species] == 1 && in_reproduction_age(agent, model)
    growth_rate = model.growthrates[agent.species]
    # W = agent.W >= 0.0 ? agent.W : 0.0
    # nchildren = rand(Poisson(growth_rate * W))
    nchildren = rand(Poisson(growth_rate))
    nchildren == 0 && return
    interaction_history = deepcopy(agent.interaction_history)
    interaction_history[1:end] .= 0
    for c in 1:nchildren
      offspring = add_agent!(agent.pos, model, agent.species, agent.biotic_phenotype, agent.abiotic_phenotype, agent.epistasisMat, agent.pleiotropyMat, agent.q, 0, agent.sex, interaction_history, model.initial_energy[agent.species], agent.W, true)
      mutate!(offspring, model)
    end
  end
end

function in_reproduction_age(agent, model)
  agent.age >= model.repro_start[agent.species] && agent.age <= model.repro_end[agent.species]
end