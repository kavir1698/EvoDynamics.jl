@testset "Mutation" begin
  dd = EvoDynamics.load_parameters(param_file)
  model = EvoDynamics.model_initiation(dd)
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

end

# #TODO
# @testset "Migration" begin

# end

# #TODO
# @testset "Reproduction" begin

# end
