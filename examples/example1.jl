
# # Simple Wright-Fisher

# We can create and run simple Wright-Fisher simulations with EvoDynamics.jl. To that end, we define a single haploid species, in an unstructured space, with two single genes affecting biotic and abiotic traits, respectively. 

using EvoDynamics

# A simple one-species model with no spatial structure. Model parameters are in a YAML file as follows:

# ```yml
# - species:
#   - id: 1
#     number of genes: 2
#     number of phenotypes: 2
#     abiotic phenotypes: [1]
#     biotic phenotypes: [2]
#     migration phenotype: 0
#     migration threshold: 3.5
#     vision radius: 0
#     check fraction: 0
#     ploidy: 1
#     epistasis matrix: [1.0, 0.0, 0.0, 1.0]
#     pleiotropy matrix: [1, 0, 0, 1]
#     growth rate: 1.0
#     expression array: [0.28, 0.46]
#     selection coefficient: 0.5
#     mutation probabilities: [0.9, 0.0, 0.0]
#     mutation magnitudes: [0.05, 0.0, 0.0]
#     N: [100]
#     environmental noise: 0.01
#     optimal phenotype values:
#       - [1.76]
#     optimal phenotypes: [1, 1, 1, 1, 1, 1]
#     age: 2
#     recombination: 0
#     initial energy: 0 

# - model:
#   generations: 5
#   space: [1,1]
#   metric: chebyshev
#   periodic: false 
#   resources: [200]
#   interactions: [-0.1] 
#   food sources: [1.0]
#   seed: Null
# ```

param_file = "../../examples/paramfile1.yml" #hide
agentdata, modeldata, model = runmodel(param_file);

modeldata