Here, I am simulating the modified version of well studied **Moran process (https://en.wikipedia.org/wiki/Moran_process).**

# Model: #

Consider a haploid population of size N, which is composed of two alleles A and a present in 1 and N - 1 copies respectively. The fitnesses of A and a are 1 + s and 1 respectively, where s > 0. 

At each discrete time step, one individual is chosen proportional to its **fitness** from the population.
It produces U number of progenies, drawn from the offspring distribution P_{U}. In order to keep the population size constant in every time step, U individuals are selected to die **randomly** from other N − 1 individuals.

P_{U} = a power law function.

 *In the Moran process, U can take only one value i.e., 2, while in this model, U is a random variable that can take any value ranging from 2 to N.*

## Aim: to calculate the probability with which the allele A will take over the population, i.e., the fixation probability of A. ##

*For more details about the model and to get an insight into the results with explanations, **see results/Results.ipynb present in this repository.***
