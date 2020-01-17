@testset "Mutation" begin
  Random.seed!(2)

  model = EvoDynamics.model_initiation(;parameters...)
  ag1 = model.agents[1]
  ag2 = model.agents[2]

  ag1.A[1] = 0.23423423
  @test ag2.A[1] != ag1.A[1]
  ag1.W = 0.123
  ag2.W = 0.321
  @test ag1.W != ag2.W
  ag1.species = 1
  ag2.species = 2
  @test ag1.species != ag2.species

  model = EvoDynamics.model_initiation(;parameters...)
  ag1 = deepcopy(model.agents[1])
  ag2 = deepcopy(model.agents[2])
  EvoDynamics.mutation!(model.agents[1], model)
  EvoDynamics.update_fitness!(model.agents[1], model)

  @test ag1.A != model.agents[1].A
  @test ag2.A == model.agents[2].A
end

@testset "Mutation and selection" begin
  Random.seed!(2)
  model = EvoDynamics.model_initiation(;parameters...)
  # We ran rounds of mutation and selection to make sure selectio doesn't affect
  # mutation
  for i in 1:00
    EvoDynamics.mutation!(model)
    EvoDynamics.selection!(model)
  end

  ag1 = deepcopy(model.agents[1])
  ag2 = deepcopy(model.agents[2])
  for i in 1:10
    EvoDynamics.mutation!(model.agents[1], model)
    EvoDynamics.update_fitness!(model.agents[1], model)
  end

  @test ag1.A != model.agents[1].A
  @test ag2.A == model.agents[2].A
end

@testset "lotkaVoltera" begin
  parameters2 = deepcopy(parameters)
  parameters2[:N] =  Dict(1 => (100, 200), 2 => (200, 100))
  parameters2[:K] = Dict(1 => [1000, 1000], 2 => [500, 500], 3 => [1000, 1000], 4 => [1000, 1000])
  parameters2[:R] = (0.1, 0.1)
  parameters2[:C] = reshape([0.1 for i in 1:4], 2, 2)
  model = EvoDynamics.model_initiation(;parameters2...)

  a11 = EvoDynamics.lotkaVoltera(model, 1, 1)
  a21 = EvoDynamics.lotkaVoltera(model, 2, 1)
  a12 = EvoDynamics.lotkaVoltera(model, 1, 2)
  a22 = EvoDynamics.lotkaVoltera(model, 2, 2)

  @test a11/100 > a21/200
  @test a12/200 < a22/100
end