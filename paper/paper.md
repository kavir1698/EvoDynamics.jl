---
title: 'EvoDynamics.jl: a framework for modeling eco-evolutionary dynamics'
tags:
  - Julia
  - evolutionary biology
  - ecology
  - agent-based modeling
  - genotype-phenotype map
  - genetic architecture
authors:
  - name: Ali R. Vahdati
    orcid: 0000-0003-0895-1495
    equal-contrib: true
    affiliation: 1
  - name: Carlos J. Melián
    orcid: 0000-0003-3974-6515
    equal-contrib: true
    affiliation: "2, 3"
affiliations:
 - name: University of Zurich, Zurich, Switzerland
   index: 1
 - name: Department of Fish Ecology and Evolution, Eawag, Centre of Ecology, Evolution and Biogeochemistry, Switzerland
   index: 2
 - name: Institute of Ecology and Evolution, Aquatic Ecology, University of Bern, Baltzerstrasse 6, CH-3012, Bern, Switzerland.
   index: 3
date: 10 September 2022
bibliography: paper.bib
---

# Summary

Genotype-phenotype maps are usually complex. Many studies have shown the mechanisms driving genome and phenome evolution mostly in isolation of one another. Yet, simulations explicitly taking into account the complex genotype-phenotype architecture are lacking. `EvoDynamics.jl` aims to connect genomes and phenomes in an easy building-block way to allow exploraing many architectures connecting the two.

Genotype-phenotype mapping usually occurs within a species. `EvoDynamics.jl` second aim is to tackle the coupling not only within a species but also between species. This is mostly because the connection between the genotype-phenotype map to biodiversity patterns is mostly unknown. This connection is particularly important for understanding the response of genetic-phenotype architectures at the biodiversity scale to global change. 

Coupling genomes to phenomes is full of challenges. This is because evolutionary dynamics, i.e., changes in allele frequencies, selective pressure on phenotypes, multi-level and antagonistic selection, among others, are functions of innumerable evolutionary and ecological processes across levels, time and space, such as interactions among genes, genotypes and phenotypes, phenotypes and abiotic environment, and complex species interactions [@loeuilleInfluence2010;@schoenerNewest2011;@ellegrenDeterminants2016].

Evolutionary modeling have mostly focused on genes and genomic processes from a variety of angles. Ecological modeling, on the other side, have mostly focused on species level interactions with other species and the environment. EvoDynamics.jl attempts to connect the two by integrating genomic and phenotypic processes within and between species. Such a connection is not easy to manage analytically. To address this complexity, in-silico simulations can provide new ideas and synthesis for a wide range of researchers and interdisciplinary teams. 

`EvoDynamics.jl` can enrich the range of applications and questions previously addressed. For example, the relationship between genetic and species diversity could not be explored accounting explicitly by the genetic-phenotypic architecture within and between each species. By accounting for the architecture within and between species, we could now explore such a connection. Similar problems like the causes of gradients in genetic diversity, the factors affecting population structure, and the effect of human actions on species and genetic variation [@leighOpportunities2021] can also be explored using `EvoDynamics.jl`. 

# Statement of need

`EvoDynamics.jl` is an open source framework for studying the link between ecological and evolutionary systems, i.e. eco-evolutionary interactions [@postEcoevolutionary2009].
The package aims to 1) make it possible for non-programmers to connect micro to macro in evolution and ecology, 2) through an agent-based model with explicit genotype-phenotype architecture, create a continuum of biological levels to understand their coupling and robustness and 3) provide a tool to ease the reproducibility, replication and extension of the results [@alstonBeginner2021]. This framework will not only make building complex simulations accessible to a large fraction of biologists who are not programmers, but also promote best practices in building models to make it accessible to other researchers.

There is a large number of frameworks for studying genetic diversity. Currently, 236 genetic simulators are listed on the [National Cancer Institute's website](https://surveillance.cancer.gov/genetic-simulation-resources/packages/). These tools can be divided into two classes of backward-time (also knows as coalescent-based) and forward-time simulators (e.g. @guillaumeNemo2006;@schiffersALADYN2014;@hallerSLiM2017;@curratSPLATCHE32019;@zhangAdmixSim2021;@bocediRangeShifter2021). Forward-time model are more powerful in handling complex evolutionary scenarios.

`EvoDynamics.jl` allows building agent-based models with complex genotype-phenotype architectures throughout providing epistasis and pleiotropy matrices as well as gene expression. In addition, the agents can be connected between different types, i.e., species, producing arbitrarily complex ecological networks, allowing any kind of species-species interactions, i.e., mutualism, commensalism, parasitism, competition. To our knowledge, this is the only framework that connects genetic architecture to such level of species interactions.

Being written in the Julia language, it is both high-performance and easily accessible to users to investigate the implementation of the code and to modify it, if needed. This is in contrast with most of the packages written in a low-level language (e.g. C++) and have an interface in an high-level language (e.g. @guillaumeNemo2006;@schiffersALADYN2014;@curratSPLATCHE32019;@bocediRangeShifter2021). Additionally, `EvoDynamics.jl` is built on top of `Agent.jl` package [@datserisAgents2022], which gives users great flexibility in defining any kind of data they want to collect. `EvoDynamics.jl` is being used in a couple of studies to understand the interplay between the genotype and phenotypic architectures to understand biodiversity patterns and species coexistence.



# Acknowledgements


# References

