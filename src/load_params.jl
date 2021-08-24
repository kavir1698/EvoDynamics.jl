function load_parameters(yaml_file::String)
  d = YAML.load_file(yaml_file)
  
  check_yml_params(d)
  check_param_shapes(d)
  reformat_params!(d)

  outd = Dict(:species => d["species"], :model => d["model"])

  return outd
end
