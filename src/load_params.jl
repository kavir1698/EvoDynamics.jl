function load_parameters(yaml_file::String)
  d = YAML.load_file(yaml_file)
  
  species_index, model_index = species_model_indices(d)

  check_yml_params(d, species_index, model_index)
  check_param_shapes(d, species_index, model_index)
  reformat_params!(d, species_index, model_index)

  outd = Dict(:species => d[species_index]["species"], :model => d[model_index])

  return outd
end

function species_model_indices(d)
  species_index = 0
  model_index = 0
  if haskey(d[1], "species") && haskey(d[2], "model")
    species_index = 1
    model_index = 2
  elseif haskey(d[2], "species") && haskey(d[1], "model")
    species_index = 2
    model_index = 1
  else
    @error "Top level species and/or model params are missing in the YAML file"
  end
  return species_index, model_index
end
