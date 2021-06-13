module EvoDynamics

using Agents
using Random
using Distributions
using StatsBase
using LinearAlgebra
using StaticArrays
using YAML

include("simulation.jl")
include("api.jl")
include("load_params.jl")
include("check_params.jl")

end # module
