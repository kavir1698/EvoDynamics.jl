function check_yml_params(d, species_index, model_index)
  species_keys = ["migration threshold", "number of genes", "number of phenotypes", "abiotic phenotypes", "biotic phenotypes", "migration phenotype", "migration threshold", "ploidy", "epistasis matrix", "pleiotropy matrix", "growth rate", "expression array", "selection coefficient", "mutation probabilities", "mutation magnitudes", "N", "vision radius", "check fraction", "environmental noise", "optimal phenotype values", "optimal phenotypes", "age", "recombination", "initial energy", "bottleneck function", "bottleneck times"]
  model_keys = ["generations", "space", "metric", "periodic", "resources", "interactions", "food sources", "seed"]
  for sp in d[species_index]["species"]
    @assert haskey(sp, "id") "Species ID field missing."
    for kk in species_keys
      @assert haskey(sp, kk) "Species $(sp["id"]) does not have field \"kk\"."
    end
  end
  for kk in model_keys 
    @assert haskey(d[model_index], kk) "Model parameter $kk is missing."
  end
end

function check_param_shapes(d, species_index, model_index)
  nspecies = length(d[species_index]["species"])

  for species in 1:nspecies
    spacesize = prod(d[model_index]["space"])
    dd = d[species_index]["species"][species]
    # Size of optimal phenotypes
    for opt in dd["optimal phenotype values"]
      @assert length(opt) == length(dd["abiotic phenotypes"]) * spacesize "Not as many optimal phenotypes as abiotic phenotypes. Species: $species"
      # for mat in opt
      #   @assert length(mat) == spacesize "Optimal phenotpyes per trait is not the same size as the number of sites"
      # end
    end
    # Biotic and abiotic phenotypes are Array
    @assert typeof(dd["abiotic phenotypes"]) <: AbstractArray "Abiotic phenotypes should be array"
    @assert typeof(dd["biotic phenotypes"]) <: AbstractArray "Biotic phenotypes should be array"
    # Ploidy
    @assert dd["ploidy"] < 3  "Ploidy more than 2 is not implemented"
    # epistasis matrix
    # @assert size(dd["epistasis matrix"], 2) % dd["ploidy"] == 0 "Number of columns in epistasisMat are not correct. They should a factor of ploidy"
    @assert length(dd["epistasis matrix"])  == (dd["ploidy"] * dd["number of genes"])^2 "epistasisMat does not have correct number of elements."
    # age is integer
    @assert typeof(dd["age"]) <: Int "Age of species $species should be integer."
    # id is integer
    @assert typeof(dd["id"]) <: Int "Age of species $species should be integer."
    # Check fraction is float
    @assert typeof(dd["check fraction"]) <: Real "Check fraction of species $species should be number."
    # recombination is Bool
    @assert typeof(dd["recombination"]) <: Real "recombination of species $species should be type numeric"
    @assert typeof(dd["initial energy"]) <: Real "Initial energy should be a number"
    # bottleneck should be nothing or array/string
    @assert typeof(dd["bottleneck function"]) <: AbstractString || isnothing(dd["bottleneck function"])
    if !isnothing(dd["bottleneck function"])
      @assert isfile(dd["bottleneck function"]) "bottleneck function points to a nonexistant file"
    end
    @assert typeof(dd["bottleneck times"]) <: AbstractArray{Int} || isnothing(dd["bottleneck times"])
  end
  # Space should be 2D
  if !isnothing(d[model_index]["space"]) && d[model_index]["space"] != "nothing"
    @assert length(d[model_index]["space"]) == 2 "Space should be only 2D"
  end
  # Resources
  @assert typeof(d[model_index]["resources"]) <: AbstractArray{Int} "Resources should be array of Integers"
  @assert length(d[model_index]["resources"]) == prod(d[model_index]["space"]) "Resources should be as long as number of sites"
  # # Area
  # @assert typeof(d[model_index]["area"]) <: AbstractArray "Area should be array"
  # @assert eltype(d[model_index]["area"]) <: Real "Area elements should be numbers"
  # @assert length(d[model_index]["area"]) == prod(d[model_index]["space"]) "Area should be as long as number of sites"
  # Interactions
  @assert typeof(d[model_index]["interactions"]) <: AbstractArray "Interactions should be an array"
  @assert eltype(d[model_index]["interactions"]) <: AbstractFloat "Elements of interactions should be floating numbers"
  @assert length(d[model_index]["interactions"]) == nspecies*nspecies "Interactions should be as many as number of species times number of species"
  # Food sources
  @assert typeof(d[model_index]["food sources"]) <: AbstractArray "Food sources should be array"
  @assert eltype(d[model_index]["food sources"]) <: AbstractFloat "Elements of food sources should be floating numbers"
  @assert length(d[model_index]["food sources"]) == nspecies*nspecies "Food sources should be as many as number of species times number of species"
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
    # Check fraction to float
    d[species_index]["species"][species]["check fraction"] = Float64(d[species_index]["species"][species]["check fraction"])
    # include bottleneck function
    if !isnothing(d[species_index]["species"][species]["bottleneck function"])
      include(d[species_index]["species"][species]["bottleneck function"])
      # reshape bottleneck times
      d[species_index]["species"][species]["bottleneck times"] = Bool.(reshape(d[species_index]["species"][species]["bottleneck times"], prod(d[model_index]["space"]), d[model_index]["generations"]))
    else
      d[species_index]["species"][species]["bottleneck times"] = Bool.(zeros(Int, prod(d[model_index]["space"]), d[model_index]["generations"]))
    end
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

  # reshape and convert resources
  d[model_index]["resources"] = reshape(d[model_index]["resources"], d[model_index]["space"]...)
  # # reshape and convert area
  # d[model_index]["area"] = reshape(d[model_index]["area"], d[model_index]["space"]...)
  # reshape and convert interactions
  d[model_index]["interactions"] = reshape(d[model_index]["interactions"], nspecies, nspecies)
  # reshape and convert food sources
  d[model_index]["food sources"] = reshape(d[model_index]["food sources"], nspecies, nspecies)
end
