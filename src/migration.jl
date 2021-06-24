
function migrate!(agent::Ind, model::ABM)
  if isnothing(model.migration_traits[agent.species]) || model.migration_thresholds[agent.species] > get_migration_trait(agent, model)  # if migration value is lower than the threshold
    return
  end

  sites = collect(nearby_positions(agent, model, model.vision_radius[agent.species]))
  if model.check_fraction[agent.species] == 0
    destination = rand(sites)
  else
    nsites = length(sites)
    nsites_selected = ceil(Int, model.check_fraction[agent.species] * nsites)
    sites_checked = EvoDynamics.sample(sites, nsites_selected, replace=false)
    # check sites and move to the best one
    distance, site_index = pick_site(agent, sites_checked, nsites_selected, model)
    current_place = check_site(agent, agent.pos, model)
    if distance < current_place
      return
    end
    destination = sites_checked[site_index]
  end
  # NB: add migration cost? If it survives migration, then do the migration.
  move_agent!(agent, destination, model)
  # agent's base fitness is its abiotic fitness 
  agent.W = abiotic_fitness(agent, model)
  return
end

"Agent evaluates the site and gives it a score"
function check_site(agent, site, model)
  phenotype = agent.abiotic_phenotype
  optimal = return_opt_phenotype(agent.species, model.step[1], site, model)
  abiotic_distance(phenotype, optimal, variance)
end

function pick_site(agent, sites, nsites, model)
  scores = Array{Float64, 1}(undef, nsites)
  for ii in 1:nsites
    scores[ii] = check_site(agent, sites[ii], model)
  end
  return findmin(scores)
end

function get_migration_trait(agent, model)
  mtrait = agent.pleiotropyMat[model.migration_traits[agent.species], :]
  sum(agent.epistasisMat[mtrait, :] * agent.q)
end
