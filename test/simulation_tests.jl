@testset "Mutation" begin
  Random.seed!(2)

  model = EvoDynamics.model_initiation(;parameters...)
  ag1 = model.agents[1]
  ag2 = model.agents[2]

  ag1.y[1] = 0.23423423
  @test ag2.y[1] != ag1.y[1]
  ag1.B[1] = false
  ag2.B[1] = true
  @test ag1.B[1] != ag2.B[1]
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
  EvoDynamics.mutation!(model.agents[1], model)
  EvoDynamics.update_fitness!(model.agents[1], model)

  @test ag1.W != model.agents[1].W
  @test ag1.B != model.agents[1].B
  @test ag1.y != model.agents[1].y
  @test ag2.W == model.agents[2].W
  @test ag2.B == model.agents[2].B
  @test ag2.y == model.agents[2].y
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

  @test ag1.W != model.agents[1].W
  @test ag1.B != model.agents[1].B
  @test ag1.y != model.agents[1].y
  @test ag2.W == model.agents[2].W
  @test ag2.B == model.agents[2].B
  @test ag2.y == model.agents[2].y
end