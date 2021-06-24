#TODO: more tests
@testset "Model" begin
  adata, mdata, model = runmodel(param_file)
  @test size(mdata, 1) == 15
  @test maximum(mdata.step) == 14
end