# EvoDynamics.jl Documentation

EvoDynamics.jl tries to bridge the gap in studying biological systems at small and large scales. Some studies only focus on single genes affecting single phenotypes, some studies only analyze gene interactions, some focus on populations, and some on species interactions. EvoDynamics.jl is a framework to study the effect of interactions between all these levels. It includes explicit pleiotropy, epistasis, selection acting on multiple phenotypes, different phenotypes affecting fitness at different amounts, arbitrary spatial structure, migration, and interacting species.

Figure below shows different biological levels controlled by the model.

![Fig. 1. __Model structure.__](struct.png)

See [Tutorial](@ref) for running the model, and [API](@ref) for a description of model parameters and simulation outline.
