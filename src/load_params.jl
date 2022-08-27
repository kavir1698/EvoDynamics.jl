function load_parameters(paramfile::String)
  include(paramfile)
  check_yml_params(model_parameters)
  check_param_shapes(model_parameters)
  return model_parameters
end
