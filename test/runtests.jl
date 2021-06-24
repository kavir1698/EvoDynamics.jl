using Test
using EvoDynamics

param_file = "params.yml"
include("simulation_tests.jl")
include("api_tests.jl")