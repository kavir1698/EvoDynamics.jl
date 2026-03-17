# EvoDynamics.jl

<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://kavir1698.github.io/EvoDynamics.jl/stable) -->
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://kavir1698.github.io/EvoDynamics.jl/dev)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/kavir1698/EvoDynamics.jl/CI.yml?branch=master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04775/status.svg)](https://doi.org/10.21105/joss.04775)

EvoDynamics.jl is a powerful framework designed to explore and understand the dynamics of biological systems across multiple scales. It offers a comprehensive set of tools and models to study evolutionary and ecological processes, providing a bridge between genotypes, phenotypes, and their intricate interactions with the environment and other species. By incorporating eco-evolutionary feedbacks, EvoDynamics.jl allows for the investigation of rapidly evolving populations that continuously adapt to changing conditions, particularly in antagonistic interactions, competitive scenarios, and mutualistic relationships.

See the [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://kavir1698.github.io/EvoDynamics.jl/dev) for more details.

## Installation

Install using the following command in Julia REPL:

```julia
]add EvoDynamics
```

## Contributions

Any contribution to EvoDynamics.jl is welcome in the following ways:

  * Modifying the code or documentation with a pull request.
  * Reporting bugs and suggestions in the issues section of the project's Github.

### Previewing Documentation Edits

Modifications to the documentation can be previewed by building the documentation locally, which is made possible by a script located in docs/make.jl. The Documenter package is required and can be installed by running `import Pkg; Pkg.add("Documenter")` in a REPL session. Then the documentation can be built and previewed in build/ first by running `julia docs/make.jl` from a terminal.

