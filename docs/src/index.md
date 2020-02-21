# EvoDynamics.jl Documentation

EvoDynamics.jl tries to bridge the gap in studying biological systems at small and large scales. Some studies only focus on single genes affecting single phenotypes, some studies only analyze gene interactions, some focus on populations, and some on species interactions. EvoDynamics.jl is a framework to study the effect of interactions between all these levels. It includes explicit pleiotropy, epistasis, selection acting on multiple phenotypes, different phenotypes affecting fitness at different amounts, arbitrary spatial structure, migration, and interacting species.

Figure below shows different biological levels controlled by the model.

![Fig. 1. __Model structure.__](struct.png)

We used the paper below as a starting point for this project. But EvoDynamics.jl goes way beyond that system by including multi-species interactions, spatial structure, and explicitly implementing epistasis and gene expression.

	Melo, D., & Marroig, G. (2015). Directional selection can drive the evolution of modularity in complex traits. Proceedings of the National Academy of Sciences, 112(2), 470–475. https://doi.org/10.1073/pnas.1322632112

See [Tutorial](@ref) for running the model, and [API](@ref) for a description of model parameters and simulation outline.
