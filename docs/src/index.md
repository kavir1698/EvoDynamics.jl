# EvoDynamics.jl Documentation

EvoDynamics.jl tries to bridge the gap in studying biological systems at small and large scales. Some studies only focus on single genes affecting single phenotypes, some studies only analyze gene interactions, some focus on populations, and some on species interactions. EvoDynamics.jl is a framework to study the effect of interactions between all these levels. It includes explicit pleiotropy, epistasis, selection acting on multiple phenotypes, different phenotypes affecting fitness at different amounts, arbitrary spatial structure, migration, and interacting species.

Figure below shows different biological levels controlled by the model.

![Fig. 1. __Model structure.__](struct.png)

See [Tutorial](@ref) for running the model, and [Model description](@ref) for a description of model parameters and simulation outline.

## Features

* Possibility to model complex food webs with various asymmetrical interactions.
* Individuals interact given their phenotypes.
* Connecting genome structure to phenotypes to populations.
* Space is a grid on which individuals can migrate.
* Possibility to model a complex environment with various amounts of resources at each site.
* Allows time varying optimal phenotypic value per site.
* Possibility of killing specific individuals at certain times and sites.
* Can model both haploid and diploid species.
* Includes mutation and recombination.
* Easy data collection.
* Runs replicates and collects data in a table automatically.
* Can run replicated in parallel.
