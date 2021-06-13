function check_yml_params(d, species_index, model_index)
  species_keys = ["migration threshold", "number of genes", "number of phenotypes", "abiotic phenotypes", "biotic phenotypes", "migration phenotype", "migration threshold", "ploidy", "epistasis matrix", "pleiotropy matrix", "covariance matrix", "growth rate", "expression array", "selection coefficient", "mutation probabilities", "mutation magnitudes", "N", "K", "vision radius", "check fraction", "environmental noise", "optimal phenotype values", "optimal phenotypes"]
  model_keys = ["generations", "space", "metric", "periodic", "seed"]
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

function check_param_shapes(d, species_index, model_index)
  nspecies = length(d[species_index]["species"])

  for species in 1:nspecies
    spacesize = prod(d[model_index]["space"])
    dd = d[species_index]["species"][species]
    # Size of optimal phenotypes
    for opt in dd["optimal phenotype values"]
      @assert length(opt) == length(dd["abiotic phenotypes"]) "Not as many optimal phenotypes as abiotic phenotypes. Species: $species"
      for mat in opt
        @assert length(mat) == spacesize "Optimal phenotpyes per trait is not the same size as the number of sites"
      end
    end
    # Biotic, abiotic and migration traits not overlapping
    theintersect = intersect(dd["biotic phenotypes"], dd["abiotic phenotypes"], dd["migration phenotype"])
    @assert length(theintersect) == 0 "There are similar traits identified as biotic and/or abiotic and/or migration"
    # Biotic and abiotic phenotypes are Array
    @assert typeof(dd["abiotic phenotypes"]) <: AbstractArray "abiotic phenotpes should be array"
    @assert typeof(dd["biotic phenotypes"]) <: AbstractArray "biotic phenotpes should be array"
    # Ploidy
    @assert dd["ploidy"] < 3  "Ploidy more than 2 is not implemented"
    # epistasis matrix
    @assert size(dd["epistasis matrix"], 2) % dd["ploidy"] == 0 "number of columns in epistasisMat are not correct. They should a factor of ploidy"
    # K and N
    @assert length(dd["K"]) == length(dd["N"]) == spacesize "There should be an N and K value for every site. Species: $species"
  end
  # Space should be 2D
  if !isnothing(d[model_index]["space"]) && d[model_index]["space"] != "nothing"
    @assert length(d[model_index]["space"]) == 2 "Space should be only 2D"
  end 
end

function reformat_params!(d, species_index, model_index)
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
    # 4. reshape optimal phenotype values
    spacesize = prod(d[model_index]["space"])
    optphens = d[species_index]["species"][species]["optimal phenotype values"]
    output = [Array{Float64, 2}[] for i in  1:length(optphens)]
    for optind in 1:length(optphens)
      start = -spacesize + 1
      for phen in 1:length(d[species_index]["species"][species]["abiotic phenotypes"])
        first = start + (phen * spacesize)
        last = first + spacesize - 1
        push!(output[optind], reshape(optphens[optind][first:last], Tuple(d[model_index]["space"])...))
      end
    end
    d[species_index]["species"][species]["optimal phenotype values"] = output
    # 5. Reshape covariance matrix
    nphenotypes = d[species_index]["species"][species]["number of phenotypes"]
    d[species_index]["species"][species]["covariance matrix"] = reshape(d[species_index]["species"][species]["covariance matrix"], nphenotypes, nphenotypes)
  end
  
  # Metric to Symbol
  d[model_index]["metric"] = Symbol(d[model_index]["metric"])

  # "nothing" to nothing for seed and space
  if d[model_index]["seed"] == "nothing"
    d[model_index]["seed"] = nothing
  end
  if d[model_index]["space"] == "nothing"
    d[model_index]["space"] = nothing
  end

  # space to tuple
  if !isnothing(d[model_index]["space"])
    d[model_index]["space"] = Tuple(d[model_index]["space"])
  end
end
