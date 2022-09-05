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
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Carlos J. Melian
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
affiliations:
 - name: University of Zurich, Zurich, Switzerland
   index: 1
 - name: Institution Name, Country
   index: 2
date: 5 September 2022p
bibliography: paper.bib
---

# Summary

Evolutionary dynamics (changes in allele frequencies, selective pressure on phenotypes, multi-level and antagonistic selection) are functions of innumerable evolutionary and ecological processes across time and space, such as interactions among genes, genotypes and phenotypes, phenotypes and abiotic environment, and complex species interactions [loeuilleInfluence2010;@schoenerNewest2011;@ellegrenDeterminants2016].
Traditionally, evolutionary theory have focused on one or a few levels of biological organization, with limited entities (e.g.genes, species) because such complex interactions are not easy to manage analytically.
To address the complexity, in silico simulations have been widely used.
These tools have been used in a wide range of applications, to answer questions such as the relationship between genetic diversity and species richness, the causes of gradients in genetic diversity, factor affecting population structure, and the effect of human actions on species and genetic variation [@leighOpportunities2021].
A problem with such simulations is reproducibility, replication and extension of the results [@alstonBeginner2021].
This is where a modeling framework becomes valuable. It will not only building complex simulations accessible to a large fraction of biologists who are not programmers, but also promote best practices in building models to make it accessible to other researchers. 

# Statement of need

`EvoDynamics.jl` is an open source framework for studying the link between ecological and evolutionary systems, i.e. eco-evolutionary interactions [@postEcoevolutionary2009].
The package aims to 1) make it possible for non-programmers to connect micro to macro in evolution and ecology, 2) through an agent-based model with explicit genetic architecture, create a continuum of biological levels of organization to understand their robustness and 3) provide a tool to ease the reproducibility and replication of results.
There is a large of number frameworks for studying genetic diversity. Currently, 236 genetic simulators are listed on the [National Cancer Institute's website](https://surveillance.cancer.gov/genetic-simulation-resources/packages/). These tools can be divided into two classes of backward-time (also knows as coalescent-based) and forward-time simulators (e.g. @guillaumeNemo2006;@schiffersALADYN2014;@hallerSLiM2017;@curratSPLATCHE32019;@zhangAdmixSim2021;@bocediRangeShifter2021).
Forward-time model are more powerful in handling complex evolutionary scenarios.

`EvoDynamics.jl` allows building agent-based models with complex genetic architectures through providing epistasis and pleiotropy matrices, and connecting that to arbitrary complex food webs, allowing any kind of species-species interactions (mutualism, commensalism, parasitism, competition). To our knowledge, this is the only framework that connects genetic architecture to such level of species interactions.
Being written in the Julia language, it is both performant and easily accessible to users to investigate the implementation of the code and to modify it, if needed. This is in contrast with most of the package that are written in a low-level language (e.g. C++) and have an interface in an high-level language (e.g. @guillaumeNemo2006;@schiffersALADYN2014;@curratSPLATCHE32019;@bocediRangeShifter2021).
Additionally, `EvoDynamics.jl` is built on top of `Agent.jl` package [@datserisAgents2022], which gives users great flexibility in defining any kind of data they want to collect.

`EvoDynamics.jl` is being used in a couple of studies to understand the effect of genetic architecture on species coexistence.

# Acknowledgements


# References

