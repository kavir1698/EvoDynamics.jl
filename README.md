# EvoDynamics.jl

<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://kavir1698.github.io/EvoDynamics.jl/stable) -->
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://kavir1698.github.io/EvoDynamics.jl/dev)
![GitHub Workflow Status](https://img.shields.io/github/workflow/status/kavir1698/EvoDynamics.jl/CI)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04775/status.svg)](https://doi.org/10.21105/joss.04775)

EvoDynamics.jl is a package for simulating the evolutionary dynamics of species in a spatially structured environment. The package allows users to define species-specific and model-specific parameters to study various aspects of population genetics, ecology, and evolution. It is a fully agent-based model, allowing individual organisms to interact with their environment and other organisms in various ways.

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

