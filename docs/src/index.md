# EvoDynamics.jl Documentation

For usage, see Tutorial.

## Theory

Effect of trait arquitecture on biodiversity.

Coupling pleiotropy and covariance matrices to explore their effect on biodiversity.

This paper describes the theoretical foundation of this project.

	Melo, D., & Marroig, G. (2015). Directional selection can drive the evolution of modularity in complex traits. Proceedings of the National Academy of Sciences, 112(2), 470–475. https://doi.org/10.1073/pnas.1322632112


## Methods

The model parameters and their role in the simulations:

* m loci.
* p traits.
* y<sub>i</sub> with i from 1 to 2m. It is representative of the m loci for a diploid individual.
* x<sub>i</sub> with i from 1 to p. It represents the additive effects of y<sub>i</sub>.
* x = B'y. If B<sub>ij</sub>=1, locus y<sub>i</sub> influences additive effect x<sub>j</sub>.
* z<sub>i</sub> = By<sub>i</sub> + ε = x<sub>i</sub> + ε. This is the phenotypic variance of individual i. ε is environmental deviation taken from a normal distribution with a uniform variance V<sub>e</sub>=0.8.
* There are two mutation rates:
  1. μ per locus per generation where y<sub>i</sub> changes from a normal distribution with mean 0 and σ=0.02.
  2. μ<sub>B</sub> per generation per entry to flip entry values of the B matrix of each individual. Each position of the B matrix is independent.
* W(z): selection = e<sup>(-1/2 ((z-θ)<sup>T</sup> ω<sup>-1</sup> (z-θ)))</sup>
* θ: multivariate fitness optimum. The rate of change of θ represents different strengths of directional selection.
* ω: covariance matrix of the selection surface.
* Reproduction:
  1. Randomly select parents with replacement proportional to their fitness.
	2. One offspring per couple.
	3. Offspring is created from gametes that include one allele from each loci and the corresponding row of the B matrix.
	4. N<sub>e</sub> constant at 5k.
* Initial parameter values from Jones et al. (15-17)

<!-- 
## Problems and solutions

### Problem 

The mean fitness of populations either goes to zero or becomes too large (infinite). Here are some possible solutions:

1. Limit the upper bound of the fitness of each individual.
2. Have different fitness functions for when there are negative values in the ω matrices and when they are all positive
3. Limit the θ optimal values within y ranges.
4. Do no use inverse of the covariance matrix. Use a matrix that always has ones on the diagonal and non-positives on the off-diagonals.

### Potential solution

1. Transform each distribution T to N(0,1)
2. Set maximum total distance to optimum Z*std(T) where Z = diagonal values cov matrix (Z=10)
3. Set fitness surface covariance mat Diag == 10; random graph from U(0,1) U(0,-1); modular U(0,1) U(0,-1)
4. Normalize fitness values sum(sum(w)) to have a max of 1 -->
