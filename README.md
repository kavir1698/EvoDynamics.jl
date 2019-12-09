# EvoDynamics

<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://kavir1698.github.io/EvoDynamics.jl/stable) -->
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://kavir1698.github.io/EvoDynamics.jl/dev)
[![Build Status](https://travis-ci.org/kavir1698/EvoDynamics.jl.svg?branch=master)](https://travis-ci.org/kavir1698/EvoDynamics.jl)


Evolutionary dynamics on multi-trait networks.

## Installation

Install using the following command inside Julia:

```julia
]add https://github.com/kavir1698/EvoDynamics.jl.git
```

It is compatible with Julia 1.3+.


## Contributions

Any contribution to EvoDynamics.jl is welcome in the following ways:

  * Modifying the code or documentation with a pull request.
  * Reporting bugs and suggestions in the issues section of the project's Github.

### Previewing Documentation Edits

Modifications to the documentation can be previewed by building the documentation locally, which is made possible by a script located in docs/make.jl. The Documenter package is required and can be installed by running `import Pkg; Pkg.add("Documenter")` in a REPL session. Then the documentation can be built and previewed in build/ first by running `julia docs/make.jl` from a terminal.

