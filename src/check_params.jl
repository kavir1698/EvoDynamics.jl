function check_yml_params(d)
  species_keys = ["migration threshold", "number of genes", "number of phenotypes", "abiotic phenotypes", "biotic phenotypes", "migration phenotype", "migration threshold", "ploidy", "epistasis matrix", "pleiotropy matrix", "growth rate", "expression array", "selection coefficient", "mutation probabilities", "mutation magnitudes", "N", "vision radius", "check fraction", "environmental noise", "optimal phenotypes", "age", "recombination", "initial energy", "bottleneck function", "reproduction start age", "reproduction end age", "functions file"]
  model_keys = ["generations", "space", "metric", "periodic", "resources", "interactions", "food sources", "seed"]
  for (k, sp) in d["species"]
    @assert haskey(sp, "name") "Species name field missing."
    for kk in species_keys
      @assert haskey(sp, kk) "Species $(sp["name"]) does not have field $(kk)."
    end
  end
  for kk in model_keys 
    @assert haskey(d["model"], kk) "Model parameter $kk is missing."
  end
end

function check_param_shapes(d)
  nspecies = length(d["species"])

  for species in 1:nspecies
    spacesize = prod(d["model"]["space"])
    dd = d["species"][species]
    # Load the function file
    @assert isfile(dd["functions file"]) "_functions file_ is nonexistant"
    include(dd["functions file"])
    # Size of optimal phenotypes
    try
      optphenotype_function = eval(Symbol(dd["optimal phenotypes"]))
    catch
      @error "Unavailable function for _optimal phenotypes_ in species $species"
    end
    # for opt in dd["optimal phenotype values"]
    #   @assert length(opt) == length(dd["abiotic phenotypes"]) * spacesize "Not as many optimal phenotypes as abiotic phenotypes. Species: $species"
    # end
    # Biotic and abiotic phenotypes are Array
    @assert typeof(dd["abiotic phenotypes"]) <: AbstractArray "Abiotic phenotypes should be array"
    @assert typeof(dd["biotic phenotypes"]) <: AbstractArray "Biotic phenotypes should be array"
    # Ploidy
    @assert dd["ploidy"] < 3 "Ploidy more than 2 is not implemented"
    # epistasis matrix
    # @assert size(dd["epistasis matrix"], 2) % dd["ploidy"] == 0 "Number of columns in epistasisMat are not correct. They should a factor of ploidy"
    @assert length(dd["epistasis matrix"]) == (dd["ploidy"] * dd["number of genes"])^2 "epistasis matrix does not have correct number of elements in species $species."
    # name is string
    @assert typeof(dd["name"]) <: AbstractString "name should be a string"
    # age is integer
    @assert typeof(dd["age"]) <: Int "Age of species $species should be integer."
    @assert typeof(dd["reproduction start age"]) <: Int "reproduction start age of species $species should be integer."
    @assert typeof(dd["reproduction end age"]) <: Int "reproduction end age of species $species should be integer."
    # Check fraction is float
    @assert typeof(dd["check fraction"]) <: Real "Check fraction of species $species should be number."
    # recombination is Bool
    @assert typeof(dd["recombination"]) <: Real "recombination of species $species should be type numeric"
    @assert typeof(dd["initial energy"]) <: Real "Initial energy should be a number"
    # bottleneck should be nothing or array/string
    @assert typeof(dd["bottleneck function"]) <: AbstractString || isnothing(dd["bottleneck function"])
    if !isnothing(dd["bottleneck function"])
      try
        bottleneck_function = eval(Symbol(dd["bottleneck function"]))
      catch
        @error "Unavailable function for _bottleneck_function_ in species $species"
      end
    end
  end
  # Space should be 2D
  if !isnothing(d["model"]["space"]) && d["model"]["space"] != "nothing"
    @assert length(d["model"]["space"]) == 2 "Space should be only 2D"
  end
  # Resources
  @assert typeof(d["model"]["resources"]) <: AbstractArray{Int} "Resources should be array of Integers"
  @assert length(d["model"]["resources"]) == prod(d["model"]["space"]) "Resources should be as long as number of sites"
  # # Area
  # @assert typeof(d["model"]["area"]) <: AbstractArray "Area should be array"
  # @assert eltype(d["model"]["area"]) <: Real "Area elements should be numbers"
  # @assert length(d["model"]["area"]) == prod(d["model"]["space"]) "Area should be as long as number of sites"
  # Interactions
  @assert typeof(d["model"]["interactions"]) <: AbstractArray "Interactions should be an array"
  @assert eltype(d["model"]["interactions"]) <: AbstractFloat "Elements of interactions should be floating numbers"
  @assert length(d["model"]["interactions"]) == nspecies*nspecies "Interactions should be as many as number of species times number of species"
  # Food sources
  @assert typeof(d["model"]["food sources"]) <: AbstractArray "Food sources should be array"
  @assert eltype(d["model"]["food sources"]) <: AbstractFloat "Elements of food sources should be floating numbers"
  @assert length(d["model"]["food sources"]) == nspecies*nspecies "Food sources should be as many as number of species times number of species"
end

function reformat_params!(d)
  nspecies = length(d["species"])

  # Fix formats and shapes
  for species in 1:nspecies
    include(d["species"][species]["functions file"])
    # 1. Reshape and convert pleiotropy matrix
    pmat = d["species"][species]["pleiotropy matrix"]
    ntraits = d["species"][species]["number of phenotypes"]
    ngenes = d["species"][species]["number of genes"] * d["species"][species]["ploidy"]
    d["species"][species]["pleiotropy matrix"] = Bool.(reshape(pmat, ntraits, ngenes))
    # 2. Reshape epistasis matrix
    d["species"][species]["epistasis matrix"] = reshape(d["species"][species]["epistasis matrix"], ngenes, ngenes)
    # 3. add key ngenes
    d["species"][species]["ngenes"] = ngenes

    d["species"][species]["optimal phenotypes"] = eval(Symbol(d["species"][species]["optimal phenotypes"]))
    # Check fraction to float
    d["species"][species]["check fraction"] = Float64(d["species"][species]["check fraction"])
    # include bottleneck function
    if !isnothing(d["species"][species]["bottleneck function"])
      bottleneck_function = eval(Symbol(d["species"][species]["bottleneck function"]))
      # reshape bottleneck times
      d["species"][species]["bottleneck function"] = bottleneck_function
    else
      d["species"][species]["bottleneck times"] = Bool.(zeros(Int, prod(d["model"]["space"]), d["model"]["generations"]))
    end
  end
  
  # Metric to Symbol
  d["model"]["metric"] = Symbol(d["model"]["metric"])

  # "nothing" to nothing for seed and space
  if d["model"]["seed"] == "nothing"
    d["model"]["seed"] = nothing
  end
  if d["model"]["space"] == "nothing"
    d["model"]["space"] = nothing
  end

  # space to tuple
  if !isnothing(d["model"]["space"])
    d["model"]["space"] = Tuple(d["model"]["space"])
  end

  # reshape and convert resources
  d["model"]["resources"] = reshape(d["model"]["resources"], d["model"]["space"]...)
  # # reshape and convert area
  # d["model"]["area"] = reshape(d["model"]["area"], d["model"]["space"]...)
  # reshape and convert interactions
  d["model"]["interactions"] = reshape(d["model"]["interactions"], nspecies, nspecies)
  # reshape and convert food sources
  d["model"]["food sources"] = reshape(d["model"]["food sources"], nspecies, nspecies)
end
