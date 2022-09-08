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
  - name: Department of Anthropology, University of Zurich, Zurich, Switzerland
    index: 1
  - name: Department of Fish Ecology and Evolution, Eawag, Centre of Ecology, Evolution and Biogeochemistry, Switzerland
    index: 2
  - name: Institute of Ecology and Evolution, Aquatic Ecology, University of Bern, Baltzerstrasse 6, CH-3012, Bern, Switzerland.
    index: 3
date: 9 September 2022
bibliography: paper.bib
---

# Summary

Genotype-phenotype maps are usually complex. Studies connecting the mechanisms driving genome and phenome evolution have been done mostly in isolation. Because the rarity of studies connecting the genotype-phenotype mapping, simulations explicitly taking into account the complex genotype-phenotype architecture are incipient [@bty197]. `EvoDynamics.jl` aims to connect genomes and phenomes in an easy building-block way to allow exploring many scenarios connecting the two.

Genotype-phenotype mapping occurs within a species, but between species interactions do occur in ecosystems. Yet, the genotype-to-phenotype-to-biodiversity connection is mostly unexplored. The second aim of `EvoDynamics.jl` is to tackle the genotype-phenotype coupling not only within a species but also between species. This is particularly important for understanding the role of the genetic-phenotype architecture for biodiversity response to global change. 

Coupling genomes to phenomes is full of challenges. This is because evolutionary dynamics, i.e., changes in allele frequencies, selective pressure on phenotypes, multi-level and antagonistic selection, among others, are functions of innumerable evolutionary and ecological processes across levels, time and space, such as interactions among genes, genotypes and phenotypes, phenotypes and abiotic environment, and complex species interactions [@loeuilleInfluence2010;@schoenerNewest2011;@ellegrenDeterminants2016;@melianetal2018].

Evolutionary modeling have mostly focused on genes and genomic processes from a variety of angles. Ecological modeling, on the other side, have mostly focused on species level interactions with other species and the environment. `EvoDynamics.jl` attempts to connect the two by integrating genomic and phenotypic processes within and between species. Such a connection is not easy to manage analytically. To address this complexity, in-silico simulations can provide new ideas and synthesis for a wide range of researchers and interdisciplinary teams aiming to join complex processes occurring in the genotype-phenotype interface to biodiversity patterns. 

`EvoDynamics.jl` can enrich the range of applications and questions previously addressed in eco-evolutionary dynamics. For example, the relationship between genetic and species diversity could not be explored accounting explicitly by the genetic-phenotypic architecture within and between each species. By accounting for the architecture within and between species, we could now explore such a connection. Similar problems like the causes of gradients in genetic and species diversity, the factors affecting population structure, and the effect of human actions on species and genetic variation [@leighOpportunities2021] can also be explored using `EvoDynamics.jl`. 
 

# Statement of need

`EvoDynamics.jl` is an open source framework for studying the link between ecological and evolutionary systems, i.e. eco-evolutionary interactions [@postEcoevolutionary2009]. The summary of the aims are: 1) to make it possible for non-programmers to connect micro to macro in evolution and ecology, 2) through an agent-based model with explicit genotype-phenotype architectures, we can couple biological levels to understand its robustness to predict biodiversity patterns, and 3) to provide a tool to ease the reproducibility, replication and extension of the results [@alstonBeginner2021]. This framework will not only make building complex simulations accessible to a large fraction of biologists who are not programmers, but also promote best practices in building reproducible models to make it accessible to other researchers.

There is a large number of frameworks for studying genetic diversity. Currently, 236 genetic simulators are listed on the [National Cancer Institute's website](https://surveillance.cancer.gov/genetic-simulation-resources/packages/). These tools can be divided into two classes of backward-time (also knows as coalescent-based) and forward-time simulators (e.g. @guillaumeNemo2006;@schiffersALADYN2014;@hallerSLiM2017;@curratSPLATCHE32019;@zhangAdmixSim2021;@bocediRangeShifter2021). Forward-time model are more powerful in handling complex evolutionary scenarios. Frameworks to explore trait dynamics also exist (e.g, @JHWUENG2020100978). 

`EvoDynamics.jl` is a forward-time simulator based on agents with complex genotype-phenotype architectures that model complex genome processes like epistasis, pleiotropy, and gene expression in the context of stochastic mutation-recombination-migration dynamics. In addition, the agents can be connected between different types, i.e., species, producing arbitrarily complex ecological networks, and any type of species-species interactions, i.e., mutualism, commensalism, parasitism, competition. To our knowledge, this is the only framework that connects genotype-phenotype maps to biodiversity, while also accounting for complex species interactions.

Being written in the Julia language, it is both high-performance and easily accessible to users to investigate the implementation of the code and to modify it, if needed. This is in contrast with most of the packages that are written in a low-level language (e.g. C++) and sometimes have an interface in an high-level language (e.g. @guillaumeNemo2006;@schiffersALADYN2014;@curratSPLATCHE32019;@bocediRangeShifter2021). Additionally, `EvoDynamics.jl` is built on top of `Agent.jl` package [@datserisAgents2022], which gives users great flexibility in defining any kind of data they want to collect. `EvoDynamics.jl` is being used in a couple of studies to understand the interplay between the genotypic and phenotypic architectures to understand biodiversity patterns and species coexistence.

# Acknowledgements

CJM acknowledges the Swiss National Science Foundation grant number IZSEZ0_183490.

# References

