# EvoDynamics.jl Documentation

Welcome to the documentation of EvoDynamics.jl, a powerful framework designed to explore and understand the dynamics of biological systems across multiple scales. EvoDynamics.jl offers a comprehensive set of tools and models to study evolutionary and ecological processes, providing a bridge between genotypes, phenotypes, and their intricate interactions with the environment and other species. By incorporating eco-evolutionary feedbacks, EvoDynamics.jl allows for the investigation of rapidly evolving populations that continuously adapt to changing conditions, particularly in antagonistic interactions, competitive scenarios, and mutualistic relationships.

With EvoDynamics.jl, you can delve into various biological levels, ranging from single genes affecting individual phenotypes to population dynamics and species interactions. The framework features an explicit genotype-phenotype architecture, capturing the complexity of pleiotropy and epistasis, selection acting on multiple phenotypes, differential fitness contributions, arbitrary spatial structures, migration, and interactions among species. Whether you are interested in modeling complex food webs, understanding the emergence of cooperation, or studying the interplay between genetic variation and ecological dynamics, EvoDynamics.jl provides a versatile platform for exploring the intricacies of biological systems.

## Biological Levels Controlled by the Model

EvoDynamics.jl enables the study of interactions at different biological levels, allowing researchers to investigate various aspects of evolutionary and ecological dynamics. The package provides capabilities to model the following levels:

1. **Genes and Phenotypes**: EvoDynamics.jl allows modeling the relationships between genes and phenotypes, including the effects of single genes on single phenotypes, gene interactions, pleiotropy, and epistasis.
2. **Population Dynamics**: The framework enables modeling population dynamics, including population growth, selection acting on multiple phenotypes, differential fitness contributions, mutation, recombination, and time-variable selection strength.
3. **Species Interactions**: EvoDynamics.jl allows for modeling interactions between species, such as competition, predation, mutualism, and parasitism. It provides the flexibility to model complex food webs with various asymmetrical interactions.
4. **Spatial Structure**: The package supports a grid-based spatial structure, enabling the modeling of individual migration and spatially explicit interactions between individuals and species.

![Fig. 1. __Model Structure illustrating different biological levels controlled by EvoDynamics.jl.__](struct.png)

## Features

EvoDynamics.jl offers a wide range of features that make it a powerful tool for studying ecological and evolutionary dynamics. Some of the key features include:

* **Complex Food Web Modeling**: The package allows for the modeling of complex food webs with various asymmetrical interactions between species.
* **Phenotype-Based Interactions**: Individuals interact with each other and the environment based on their phenotypes, allowing for realistic ecological dynamics.
* **Genotype-Phenotype Mapping**: EvoDynamics.jl provides the ability to connect genome structure to phenotypes and populations, allowing for the investigation of genotype-phenotype relationships.
* **Spatial Modeling**: The package supports a grid-based spatial structure, enabling the modeling of individual migration and spatially explicit interactions.
* **Dynamic Environment**: EvoDynamics.jl allows for the modeling of a complex environment with spatio-temporally varying resources, providing a realistic setting for ecological interactions.
* **Optimal Phenotypic Values**: Researchers can incorporate time-varying optimal phenotypic values per site and per species, allowing for the study of phenotypic adaptation to changing environmental conditions.
* **Differential Fitness Contributions**: The framework supports the differential contributions of traits to fitness, allowing for the investigation of trait-dependent selection.
* **Selective Removal**: EvoDynamics.jl provides the ability to remove individuals at specific times and locations, allowing for the modeling of selective pressures and perturbations in the system.
* **Haploid and Diploid Species**: The package supports modeling both haploid and diploid species, providing flexibility in representing different reproductive modes.
* **Mutation and Recombination**: EvoDynamics.jl incorporates mutation and recombination mechanisms, allowing for the simulation of genetic variation and evolution.
* **Data Collection and Analysis**: The package offers easy data collection and analysis, allowing researchers to collect and analyze simulation results efficiently.
* **Parallel Computing**: Researchers can run replicates in parallel, leveraging the computational power of multiple processors for increased efficiency.


For detailed information on how to use EvoDynamics.jl and its various functionalities, please refer to the [Tutorial](@ref) and [Model Parameters and Simulation Outline](@ref) sections.
