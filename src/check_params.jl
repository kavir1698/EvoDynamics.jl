function check_yml_params(d, species_index, model_index)
  species_keys = ["migration threshold", "number of genes", "number of phenotypes", "abiotic phenotypes", "biotic phenotypes", "migration phenotype", "migration threshold", "ploidy", "epistasis matrix", "pleiotropy matrix", "covariance matrix", "growth rate", "expression array", "selection coefficient", "mutation probabilities", "mutation magnitudes", "N", "K", "vision radius", "check fraction", "environmental noise", "optimal phenotype values", "optimal phenotypes"]
  model_keys = ["generations", "space", "metric"]
  for sp in d[species_index]["species"]
    @assert haskey(sp, "id") "species ID field missing"
    for kk in species_keys
      @assert haskey(sp, kk) "Input error: species $(sp["id"]) does not have field \"kk\""
    end
  end
  for kk in model_keys 
    @assert haskey(d[model_index], kk) "model parameter $kk is missing"
  end
end

#TODO
function check_param_shapes()
  
end