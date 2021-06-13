function load_parameters(yaml_file::String)
  d = YAML.load_file(yaml_file)
  
  if haskey(d[1], "species") && haskey(d[2], "model")
    species_index = 1
    model_index = 2
  elseif haskey(d[2], "species") && haskey(d[1], "model")
    species_index = 2
    model_index = 1
  else
    @error "Top level species and/or model params are missing in the YAML file"
  end

  check_yml_params(d, species_index, model_index)
  
  nspecies = length(d[species_index]["species"])

  # Fix formats and shapes
  for species in 1:nspecies
    # 1. Reshape and convert pleiotropy matrix
    pmat = d[species_index]["species"][species]["pleiotropy matrix"]
    ntraits = d[species_index]["species"][species]["number of phenotypes"]
    ngenes = d[species_index]["species"][species]["number of genes"] * d[species_index]["species"][species]["ploidy"]
    d[species_index]["species"][species]["pleiotropy matrix"] = Bool.(reshape(pmat, ntraits, ngenes))
    # 2. Reshape epistasis matrix
    d[species_index]["species"][species]["epistasis matrix"] = reshape(d[species_index]["species"][species]["epistasis matrix"], ngenes, ngenes)
    # 3. add key ngenes
    d[species_index]["species"][species]["ngenes"] = ngenes
  end
  
  d[model_index]["metric"] = Symbol(d[model_index]["metric"])

  outd = Dict(:species => d[species_index]["species"], :model => d[model_index])

  return outd
end
