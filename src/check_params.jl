function check_yml_params(d)
  species_keys = [:name, :number_of_genes, :number_of_phenotypes, :abiotic_phenotypes, :biotic_phenotypes, :migration_phenotype, :migration_threshold, :ploidy, :epistasis_matrix, :pleiotropy_matrix, :growth_rate, :expression_array, :selection_coefficient, :mutation_probabilities, :mutation_magnitudes, :N, :vision_radius, :check_fraction, :environmental_noise, :optimal_phenotypes, :age, :recombination, :initial_energy, :bottleneck_function, :reproduction_start_age, :reproduction_end_age]
  model_keys = [:species, :generations, :space, :metric, :periodic, :resources, :interactions, :food_sources, :seed]
  species = d[:species]
  for sp in species
    spfields = collect(keys(sp))
    @assert in(:name, spfields) "Species struct has no field `name` "
    spname = sp[:name]
    @assert length(sp) == length(species_keys)  "Incorrect number of keys for species $spname"
    for kk in species_keys
      @assert in(kk, spfields) "Species $spname does not have field $(kk)."
    end
  end
  mfields = collect(keys(d))
  for kk in model_keys
    @assert in(kk, model_keys) "Model parameter $kk is missing."
  end
  
end

function check_param_shapes(d)
  allspecies = d[:species]
  nspecies = length(allspecies)

  for dd in allspecies
    spname = dd[:name]
    @assert typeof(dd[:name]) <: AbstractString "name should be a string"
    @assert typeof(dd[:abiotic_phenotypes]) <: AbstractArray "Abiotic phenotypes should be array"
    @assert typeof(dd[:biotic_phenotypes]) <: AbstractArray "Biotic phenotypes should be array"
    @assert eltype(dd[:pleiotropy_matrix]) <: Bool "pleiotropy matrix should contain boolean values in species $spname"
    # Ploidy
    @assert dd[:ploidy] < 3 "Ploidy more than 2 is not implemented"
    # epistasis matrix
    nloci = dd[:ploidy] * dd[:number_of_genes]
    @assert size(dd[:epistasis_matrix]) == (nloci, nloci) "epistasis matrix does not have correct size in species $spname. It should be (nloci, nloci) where nloci is ploidy * number of genes"
    @assert typeof(dd[:age]) <: Int "Age of species $spname should be integer."
    @assert typeof(dd[:reproduction_start_age]) <: Int "reproduction start age of species $spname should be integer."
    @assert typeof(dd[:reproduction_end_age]) <: Int "reproduction end age of species $spname should be integer."
    @assert typeof(dd[:check_fraction]) <: AbstractFloat "Check fraction of species $spname should be floating point number."
    # recombination is Bool
    @assert typeof(dd[:recombination]) <: Real "recombination of species $spname should be type numeric"
    @assert typeof(dd[:initial_energy]) <: Real "Initial energy should be a number"
    # bottleneck should be nothing or array/string
    @assert typeof(dd[:bottleneck_function]) <: Function || isnothing(dd[:bottleneck_function]) "bottleneck function for species $spname is not a function or `nothing`"
    @assert typeof(dd[:selection_coefficient]) <: AbstractFloat "selection coefficient of species $spname should be floating point number"
    @assert typeof(dd[:migration_threshold]) <: AbstractFloat "migration_threshold of species $spname should be floating point number"
  end
  # Space should be 2D
  if !isnothing(d[:space])
    @assert length(d[:space]) == 2 "Space should be only 2D"
  end
  @assert typeof(d[:space]) <: Tuple "Space should be either a tuple or `nothing`"
  # Resources
  
  @assert typeof(d[:resources]) <: Function "Resources function not defined"

  # Interactions
  @assert typeof(d[:interactions]) <: AbstractArray "Interactions should be an array(matrix)"
  @assert eltype(d[:interactions]) <: AbstractFloat "Elements of interactions should be floating numbers"
  @assert length(d[:interactions]) == nspecies * nspecies "Interactions should be as many as number of species times number of species"
  # Food sources
  @assert typeof(d[:food_sources]) <: AbstractArray "Food sources should be array"
  @assert eltype(d[:food_sources]) <: AbstractFloat "Elements of food sources should be floating numbers"
  @assert length(d[:food_sources]) == nspecies * nspecies "Food sources should be as many as number of species times number of species"
end
