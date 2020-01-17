var documenterSearchIndex = {"docs":
[{"location":"api/#API-1","page":"API","title":"API","text":"","category":"section"},{"location":"api/#Parameters-1","page":"API","title":"Parameters","text":"","category":"section"},{"location":"api/#","page":"API","title":"API","text":"The parameters below are required for any simulation. They should be in a dictionary object. The dictionary keys should be the parameters names as Symbol (with a colon \":\" before the name). See the Tutorial page for an example of the parameters dictionary.","category":"page"},{"location":"api/#","page":"API","title":"API","text":"m: A tuple or an array specifying the ploidy of each species. Currently only support m=1 (haploid) and m=2 (diploid).\nP: A tuple or an array specifying the number of traits p for each species.\nL: A tuple or an array specifying the number of loci l for each species.\nR: A tuple or an array specifying growth for each species. Growth rates are for a logistic growth model, where N_t+1 = N_t + rtimes Ntimes (1 - ((NK)), where N is population size, t is time, r is growth rate and K is carrying capacity. If r=0, population size remains constant.\nC: A matrix containing competition coefficients between each pair of species. A competition coefficient denotes the strength of competition exerted by an individual of species j on an individual of species i. It uses the Lotka-Voltera equation Ni_t+1 = Ni_t + rtimes Ntimes (1 - ((Ni + cNj)K) where c is competition coefficient. When competition coefficient is positive, population j competes with population i. If negative, population j helps population i to grow. And if 0, population j does not affect population i. If c_ij  0 and c_ji  0, both populations are in competition, if c_ij  0 and c_ji  0, species i is a parasite of species j. If c_ij  0 and c_ji  0, the two species have a mutualistic relationship. If c_ij  0 and c_ji = 0, they have a commensal relationship.\nA: A tuple or an array of genotype-phenotype matrices, one for each species. Each matrix shows the amount of contribution of each locus to each phenotype. It has P rows and L columns. Make sure that a row is not all zero (a trait is controlled by no locus).\nY: A tuple or an array  of selection coefficients for each species.\nT: A tuple or an array of arrays, where each inner array specifies optimal phenotypes θ for each species. Each inner array should be of length p (number of traits) of its corresponding species.\nΩ: A tuple or an array of matrices, each of which ω represents a covariance matrix of the selection surface. Each matrix is of size ptimes p.\nM: A tuple or an array of mutation rates μ for each species.\nN: A dictionary where each key is a node number and its value is a tuple for population size of each species at that node. This dictionary does not have to have a key for all the nodes, but it should have a value for all the species.\nK: A dictionary where each key is a node number and its value is tuple of carrying capacities K of the node for each species. The dictionary should have a key for all the nodes and carrying capacity for each species per node.\nmigration_rates: An array of matrices, each matrix shows of migration rates between each pair of nodes for a species. The rows and columns of the matrix are node numbers in order. If instead of a matrix, there is nothing, no migration occurs for that species.\nE: A tuple  or an array of the variance of a normal distribution ε representing environmental noise for each species.\ngenerations: number of generations to run the simulation.\nspace: Either a tuple of size 2 or 3 for a grid size or a SimpleGraph object for an arbitrary graph. If it is a tuple, a grid is built internally\nmoore: Whether nodes in the grid have 8 neighbors (Moore neighborhood). Default is false, i.e. cells only have 4 neighbors.","category":"page"},{"location":"api/#Simulation-outline-1","page":"API","title":"Simulation outline","text":"","category":"section"},{"location":"api/#","page":"API","title":"API","text":"Within each time-step, the following occurs:","category":"page"},{"location":"api/#","page":"API","title":"API","text":"Mutation\nFitness update\nMigration\nReproduction (only for diploids)\nselection","category":"page"},{"location":"api/#Mutation-1","page":"API","title":"Mutation","text":"","category":"section"},{"location":"api/#","page":"API","title":"API","text":"The genotype vector y and pleiotropy matrix B of each individual mutates.","category":"page"},{"location":"api/#","page":"API","title":"API","text":"y mutates by adding random numbers from a normal distribution with mean 0 and standard deviation from the δ to the current values of y. δ is specified by the M parameter.","category":"page"},{"location":"api/#","page":"API","title":"API","text":"B mutates by randomly switching 0s and 1s with probability given in parameter MB.","category":"page"},{"location":"api/#Fitness-update-1","page":"API","title":"Fitness update","text":"","category":"section"},{"location":"api/#","page":"API","title":"API","text":"Fitness of each individual updates after mutation. Fitness is W = exp(γ times transpose(z - θ)times inv(ω)times (z - θ)), where is the phenotype vector (z = By + μ), γ is selection coefficient, θ is optimum phenotypes vector, and ω is covariance matrix of selection surface. ","category":"page"},{"location":"api/#Migration-1","page":"API","title":"Migration","text":"","category":"section"},{"location":"api/#","page":"API","title":"API","text":"Each agent moves with probabilities given in migration_rates to another node.","category":"page"},{"location":"api/#Reproduction-1","page":"API","title":"Reproduction","text":"","category":"section"},{"location":"api/#","page":"API","title":"API","text":"When a species is diploid, they sexually reproduce. To that end, individuals of the same species in the same location are randomly paired. Each pair produces one offspring. Then the parents die.","category":"page"},{"location":"api/#","page":"API","title":"API","text":"To produce an offspring, each parent contributes to half of the offspring's genotype y and pleiotropy matrix B. The genes coming from each parent are randomly assigned.","category":"page"},{"location":"api/#Selection-1","page":"API","title":"Selection","text":"","category":"section"},{"location":"api/#","page":"API","title":"API","text":"A number of individuals n are selected for the next generation via sampling with replacement weighted by individuals' fitness values. n is calculated using the Lotka-Voltera model for population competition Ni_t+1 = Ni_t + rtimes Ntimes (1 - ((Ni + cNj)K) where N is population size, t is time, r is growth rate and K is carrying capacity, and c is competition coefficient. Briefly, each population growth with a logistic model when it is not affected by other species. Otherwise, its growth increases or decreases depending on its interactions with other species.","category":"page"},{"location":"api/#Data-collection-1","page":"API","title":"Data collection","text":"","category":"section"},{"location":"api/#","page":"API","title":"API","text":"TODO","category":"page"},{"location":"example2/#","page":"Example 2","title":"Example 2","text":"EditURL = \"https://github.com/kavir1698/EvoDynamics.jl/blob/master/examples/example2.jl\"","category":"page"},{"location":"example2/#Dummy-example-2-1","page":"Example 2","title":"Dummy example 2","text":"","category":"section"},{"location":"example2/#","page":"Example 2","title":"Example 2","text":"b = 2","category":"page"},{"location":"example2/#","page":"Example 2","title":"Example 2","text":"","category":"page"},{"location":"example2/#","page":"Example 2","title":"Example 2","text":"This page was generated using Literate.jl.","category":"page"},{"location":"example1/#","page":"Example 1","title":"Example 1","text":"EditURL = \"https://github.com/kavir1698/EvoDynamics.jl/blob/master/examples/example1.jl\"","category":"page"},{"location":"example1/#Dummy-example-1","page":"Example 1","title":"Dummy example","text":"","category":"section"},{"location":"example1/#","page":"Example 1","title":"Example 1","text":"a = 2","category":"page"},{"location":"example1/#","page":"Example 1","title":"Example 1","text":"","category":"page"},{"location":"example1/#","page":"Example 1","title":"Example 1","text":"This page was generated using Literate.jl.","category":"page"},{"location":"#EvoDynamics.jl-Documentation-1","page":"Introduction","title":"EvoDynamics.jl Documentation","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"For usage, see Tutorial.","category":"page"},{"location":"#Theory-1","page":"Introduction","title":"Theory","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"Effect of trait arquitecture on biodiversity.","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"Coupling pleiotropy and covariance matrices to explore their effect on biodiversity.","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"This paper describes the theoretical foundation of this project.","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"Melo, D., & Marroig, G. (2015). Directional selection can drive the evolution of modularity in complex traits. Proceedings of the National Academy of Sciences, 112(2), 470–475. https://doi.org/10.1073/pnas.1322632112","category":"page"},{"location":"#Methods-1","page":"Introduction","title":"Methods","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"The model parameters and their role in the simulations:","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"m loci.\np traits.\ny<sub>i</sub> with i from 1 to 2m. It is representative of the m loci for a diploid individual.\nx<sub>i</sub> with i from 1 to p. It represents the additive effects of y<sub>i</sub>.\nx = B'y. If B<sub>ij</sub>=1, locus y<sub>i</sub> influences additive effect x<sub>j</sub>.\nz<sub>i</sub> = By<sub>i</sub> + ε = x<sub>i</sub> + ε. This is the phenotypic variance of individual i. ε is environmental deviation taken from a normal distribution with a uniform variance V<sub>e</sub>=0.8.\nThere are two mutation rates:\nμ per locus per generation where y<sub>i</sub> changes from a normal distribution with mean 0 and σ=0.02.\nμ<sub>B</sub> per generation per entry to flip entry values of the B matrix of each individual. Each position of the B matrix is independent.\nW(z): selection = e<sup>(-1/2 ((z-θ)<sup>T</sup> ω<sup>-1</sup> (z-θ)))</sup>\nθ: multivariate fitness optimum. The rate of change of θ represents different strengths of directional selection.\nω: covariance matrix of the selection surface.\nReproduction:\nRandomly select parents with replacement proportional to their fitness.","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"2. One offspring per couple.\n3. Offspring is created from gametes that include one allele from each loci and the corresponding row of the B matrix.\n4. N<sub>e</sub> constant at 5k.","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"Initial parameter values from Jones et al. (15-17)","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"<!– ","category":"page"},{"location":"#Problems-and-solutions-1","page":"Introduction","title":"Problems and solutions","text":"","category":"section"},{"location":"#Problem-1","page":"Introduction","title":"Problem","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"The mean fitness of populations either goes to zero or becomes too large (infinite). Here are some possible solutions:","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"Limit the upper bound of the fitness of each individual.\nHave different fitness functions for when there are negative values in the ω matrices and when they are all positive\nLimit the θ optimal values within y ranges.\nDo no use inverse of the covariance matrix. Use a matrix that always has ones on the diagonal and non-positives on the off-diagonals.","category":"page"},{"location":"#Potential-solution-1","page":"Introduction","title":"Potential solution","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"Transform each distribution T to N(0,1)\nSet maximum total distance to optimum Z*std(T) where Z = diagonal values cov matrix (Z=10)\nSet fitness surface covariance mat Diag == 10; random graph from U(0,1) U(0,-1); modular U(0,1) U(0,-1)\nNormalize fitness values sum(sum(w)) to have a max of 1 –>","category":"page"},{"location":"tutorial/#Tutorial-1","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/#EvoDynamics.jl's-basic-usage-1","page":"Tutorial","title":"EvoDynamics.jl's basic usage","text":"","category":"section"},{"location":"tutorial/#","page":"Tutorial","title":"Tutorial","text":"First, define your model parameters. Below is a set of random parameters. See API for a description of each parameter.","category":"page"},{"location":"tutorial/#","page":"Tutorial","title":"Tutorial","text":"using Random\nimport LinearAlgebra: Symmetric\n\nP = (4, 5)\nL = (7, 8)\nm = (2, 1)\nparameters = Dict(\n  :L => L .* m,\n  :P => P,\n  :R => (0.8, 0.12),\n  :C => rand(-0.1:0.01:0.1, 2, 2),\n  :A =>  Tuple([Random.rand(i[1], i[2]) for i in zip(P, L .* m)]),\n  :Y => (-0.5, -0.5),\n  :m => m,\n  :T => Tuple([randn(Float16, n) for n in P]),\n  :Ω => Tuple([Symmetric(rand(Float16, i[1], i[2])) for i in zip(P, P)]),\n  :M => (0.02, 0.02),\n  :N => Dict(1 => (1000, 1000)),\n  :K => Dict(1 => [1000, 1000], 2 => [1000, 1000], 3 => [1000, 1000], 4 => [1000, 1000]),\n  :migration_rates => [[1.0 0.02 0.02 0.02; 0.03 1.0 0.03 0.03; 0.01 0.01 1.0 0.01; 0.01 0.01 0.01 1.0] for i in 1:2],\n  :E => (0.08, 0.08),\n  :generations => 5,\n  :space => (2,2),\n  :moore => false\n)","category":"page"},{"location":"tutorial/#","page":"Tutorial","title":"Tutorial","text":"We can the use the runmodel function to create a model from these parameters and run the simulation.","category":"page"},{"location":"tutorial/#","page":"Tutorial","title":"Tutorial","text":"EvoDynamics.runmodel","category":"page"},{"location":"tutorial/#EvoDynamics.runmodel","page":"Tutorial","title":"EvoDynamics.runmodel","text":"runmodel(parameters::Dict; kwargs)\n\nCreates and runs a model given parameters. Returns a DataFrame of collected data, which are specified by kwargs.\n\nKeywords\n\ncollect::Dict=Dict(:model => [meanfitnessper_species]) Data to be collected. By default, collects mean population fitness per species.\nwhen::AbstractArray{Int}=1:parameters[:generations] The generations from which data are collected\nreplicates::Int = 0 Number of replicates per simulation.\nparallel::Bool = false Whether to run replicates in parallel. If true, you should add processors to your julia session (e.g. by addprocs(n)) and define your parameters and EvoDynamics on all workers. To do that, add @everywhere before them. For example, @everywhere EvoDynamics.\n\n\n\n\n\n","category":"function"},{"location":"tutorial/#","page":"Tutorial","title":"Tutorial","text":"using EvoDynamics\ndata = runmodel(parameters)\ndata[1:5, :]","category":"page"}]
}
