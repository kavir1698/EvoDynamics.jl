Random.seed!(10)

@testset "Model" begin
  parameters2 = deepcopy(parameters)
  adata, mdata, model = runmodel(parameters2)
  @test size(mdata, 1) == parameters2[:generations] + 1
  @test maximum(mdata.step) == parameters2[:generations]
  
  parameters2[:space] = nothing
  parameters2[:migration_rates] = [nothing, nothing]
  adata, mdata, model = runmodel(parameters2)
  @test size(mdata, 1) == parameters2[:generations] + 1
  @test maximum(mdata.step) == parameters2[:generations]
end