@testset "Basic functions" begin
  dd = EvoDynamics.load_parameters(param_file)
  model = EvoDynamics.model_initiation(dd)
  EvoDynamics.remove_agent!(model[1], model)
  @test !haskey(model.agents, 1) == true

  EvoDynamics.consume_food!(model[2], model)
  @test model.resources == [1999 1980; 1830 1900]

  @test EvoDynamics.coord2vertex((2, 1), model) == 2
  @test EvoDynamics.coord2vertex((1, 2), model) == 3
  @test EvoDynamics.coord2vertex((2, 2), model) == 4
  @test EvoDynamics.vertex2coord(2, model) == (2, 1)
  @test EvoDynamics.vertex2coord(3, model) == (1, 2)
  @test EvoDynamics.vertex2coord(4, model) == (2, 2)

  @test model[2].W == model[3].W
  model[2].q[1] -= 0.05
  EvoDynamics.update_fitness!(model[2], model)
  @test model[2].W != model[3].W
end

@testset "Mutation" begin
  dd = EvoDynamics.load_parameters(param_file)
  model = EvoDynamics.model_initiation(dd)

  @test EvoDynamics.nagents(model) == 1100
  @test size(model.space) == (2, 2)
  @test length(EvoDynamics.ids_in_position((1, 1), model)) == EvoDynamics.nagents(model)

  ag1 = model.agents[1]
  ag2 = model.agents[2]

  ag1.epistasisMat[2] = 0.36
  @test ag2.epistasisMat[2] != ag1.epistasisMat[1]
  ag1.W = 0.123
  ag2.W = 0.321
  @test ag1.W != ag2.W
  ag1.species = 1
  ag2.species = 2
  @test ag1.species != ag2.species

  Worg = copy(ag1.W)
  ag1.q .+= 0.5
  EvoDynamics.update_fitness!(ag1, model)
  @test ag1.W != Worg

  # no leakage. using deepcopy for creating new agents
  @test model[3].q == model[4].q
  @test model[3].pleiotropyMat == model[4].pleiotropyMat
  @test model[3].epistasisMat == model[4].epistasisMat
  EvoDynamics.mutate!(model[3], model)
  @test model[3].q != model[4].q
  @test model[3].pleiotropyMat != model[4].pleiotropyMat
  @test model[3].epistasisMat != model[4].epistasisMat
  @test model[5].q == model[4].q
  @test model[5].pleiotropyMat == model[4].pleiotropyMat
  @test model[5].epistasisMat == model[4].epistasisMat
end

@testset "Seeding" begin
  dd = EvoDynamics.load_parameters(param_file)
  dd[:seed] = 1
  model = model_initiation(dd)
  EvoDynamics.step!(model, EvoDynamics.agent_step!, 1)

  dd2 = EvoDynamics.load_parameters(param_file)
  dd2[:seed] = 2
  model2 = EvoDynamics.model_initiation(dd2)
  EvoDynamics.step!(model2, EvoDynamics.agent_step!, 1)

  dd3 = EvoDynamics.load_parameters(param_file)
  dd3[:seed] = 1
  model3 = EvoDynamics.model_initiation(dd3)
  EvoDynamics.step!(model3, EvoDynamics.agent_step!, 1)

  @test collect(keys(model.agents)) != collect(keys(model2.agents))
  @test collect(keys(model.agents)) == collect(keys(model3.agents))
end

@testset "Interactions" begin
  dd = EvoDynamics.load_parameters(param_file)
  model = EvoDynamics.model_initiation(dd)

  agent = model[1]
  bph = EvoDynamics.get_biotic_phenotype(agent, model)
  aph = EvoDynamics.get_abiotic_phenotype(agent, model)

  @test length(bph) == length(dd[:species][1][:biotic_phenotypes])
  @test length(aph) == length(dd[:species][1][:abiotic_phenotypes])

  # identical phenotypes in the same species results in the same distance:
  agent2 = model[1001]
  d = EvoDynamics.phenotypic_distance(agent, agent2, model)
  dsame = EvoDynamics.phenotypic_distance(agent, agent, model)
  dsame2 = EvoDynamics.phenotypic_distance(agent2, agent2, model)
  @test dsame == dsame2 == 0

  # identical phenotypes in two different species results in non-zero distance:
  agent.biotic_phenotype .= agent2.biotic_phenotype
  d2 = EvoDynamics.phenotypic_distance(agent, agent2, model)
  @test d2 > 0.0
  @test EvoDynamics.phenotypic_distance(agent, agent2, model) == EvoDynamics.phenotypic_distance(agent2, agent, model)

  optphen = EvoDynamics.return_opt_phenotype(agent2.species, agent.pos, model)
  @test EvoDynamics.abiotic_distance(agent2.abiotic_phenotype, optphen, model.abiotic_variances[agent2.species]) > 0
  agent2.abiotic_phenotype .= optphen
  @test EvoDynamics.abiotic_distance(agent2.abiotic_phenotype, optphen, model.abiotic_variances[agent2.species]) == 0
  @test EvoDynamics.abiotic_fitness(agent2, model) == 1.0

  EvoDynamics.interact!(agent, agent2, model)
  @test agent.isalive == false

  @test length(EvoDynamics.target_species_ids(agent, model)) == 2
  @test length(EvoDynamics.target_species_ids(agent2, model)) == 2

  # reproduction
  nagentsbefore = EvoDynamics.nagents(model)
  EvoDynamics.reproduce!(agent2, model)
  EvoDynamics.reproduce!(agent2, model)
  EvoDynamics.reproduce!(agent2, model)
  nagentsafter = EvoDynamics.nagents(model)
  @test nagentsbefore < nagentsafter

  agent3 = model[10]
  agent4 = model[11]
  agent3.sex = true
  agent4.sex = false

  nagentsbefore = EvoDynamics.nagents(model)
  EvoDynamics.reproduce!(agent3, agent4, model)
  EvoDynamics.reproduce!(agent3, agent4, model)
  nagentsafter = EvoDynamics.nagents(model)
  @test nagentsbefore < nagentsafter
end
