"""
    load_parameters(paramfile::String)

Reads the parameters from a parameter file (`.jl`), checks them for correct format, and returns the `model_parameters` dictionary object. See the online docs for the list of parameters and their correct values.

The output of this function can be loaded to the `model_initiation` function to create an ABM object. A user would normally not need to do these steps directly, because they are handled by the `runmodel` function. These will be useful for debugging and running the model in a different way.
"""
function load_parameters(paramfile::String)
  include(paramfile)
  check_param_names(model_parameters)
  check_param_shapes(model_parameters)
  return model_parameters
end
