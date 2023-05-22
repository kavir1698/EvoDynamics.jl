module EvoDynamics

using Agents
using Random
using Distributions
using StatsBase
using LinearAlgebra
using StaticArrays
using YAML
using RandomNumbers

export runmodel, model_initiation, load_parameters
# export abiotic_survive!, burn_energy!, consume_food!, interact!, kill_agent!, reproduce!, migrate!, target_species_ids, eat!, abiotic_distance, abiotic_fitness, interaction_power, phenotypic_distance, get_abiotic_phenotype, get_biotic_phenotype, update_fitness!, vertex2coord, coord2vertex, mutate!, adjust_fitness!, create_gamete, crossing_overs, return_opt_phenotype
# export check_site, pick_site, get_migration_trait
# export step!

include("simulation.jl")
include("interactions.jl")
include("migration.jl")
include("reproduction.jl")
include("api.jl")
include("load_params.jl")
include("check_params.jl")

end # module
