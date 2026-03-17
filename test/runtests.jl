using Test
using EvoDynamics

param_file = "params.jl"
include("simulation_tests.jl")
include("api_tests.jl")