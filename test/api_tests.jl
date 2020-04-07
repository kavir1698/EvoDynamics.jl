Random.seed!(10)

@testset "Model" begin
  parameters2 = deepcopy(parameters)
  a, model = runmodel(parameters2)
  @test size(a, 1) == parameters2[:generations] + 1
  @test maximum(a.step) == parameters2[:generations]
  
  parameters2[:space] = nothing
  parameters2[:migration_rates] = [nothing, nothing]
  a, model = runmodel(parameters2)
  @test size(a, 1) == parameters2[:generations] + 1
  @test maximum(a.step) == parameters2[:generations]
end