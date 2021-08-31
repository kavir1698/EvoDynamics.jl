#TODO: more tests
@testset "Model" begin
  adata, mdata, models = runmodel(param_file)
  @test size(mdata, 1) == 15
  @test maximum(mdata.step) == 14
end

@testset "replicates" begin
  replicates = 3
  adata, mdata, models = runmodel(param_file, replicates = replicates)
  size(mdata)
  @test size(mdata, 1) == 15*replicates
  @test maximum(mdata.step) == 14
  @test maximum(mdata.ensemble) == 3
end