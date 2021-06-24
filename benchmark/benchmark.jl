using EvoDynamics
using BenchmarkTools

param_file = "../test/params.yml" 

a = @benchmark runmodel(param_file)
display(a)