# EvoDynamics.jl Documentation

EvoDynamics.jl tries to bridge the gap in studying biological systems at a narrow scale. Some studies only focus on single genes affecting single phenotypes, some studies only analyze gene interactions, some focus on population level. EvoDynamics.jl is a framework that make it easy to simulate a system with different levels from pleiotropy, epistasis, selection acting on multiple phenotypes, spatially structured populations, migration, and interactions between species.

We used the paper below as a starting point for this project. But EvoDynamics.jl goes way beyond that system by including multi-species interactions, spatial structure, and explicitly implementing epistasis and gene expression.

	Melo, D., & Marroig, G. (2015). Directional selection can drive the evolution of modularity in complex traits. Proceedings of the National Academy of Sciences, 112(2), 470–475. https://doi.org/10.1073/pnas.1322632112

See [Tutorial](@ref) for running the model, and [API](@ref) for a description of model parameters and simulation outline.
