# v0.17

* Update interaction rules so that the probability within species interactions are independent of the abundance of other species in a site.

# v0.16

* Users can modify the sequence of events in the simulations by providing their own stepping functions in the `runmodel` functions.
* Allows time-variable selection coefficients per species so that selection can be relaxed or increase during the simulation.
* Allow differential contribution of abiotic traits to fitness.

# v0.15

Breaking. Because of problems with using a function in the parameter file, this version uses a .jl file for parameters. Functions cannot be parameters, use arrays instead.

# v0.13.0

Breaking. Users now give functions for optimal phenotypes instead of a list of times and values. This is more flexible and concise. Same is true for bottlenecks.

# v0.11.0

Updates the internals of the API.

# v0.10.0

Includes more parameters and bug fixes.

# v0.5.0

## Breaking changes

* Moved to Agents.jl v3.0. Data collection method has changed.
* Parameter names have changed to be more verbose.

# v0.4.1

Changelog is kept with respect to version 0.4.1.
