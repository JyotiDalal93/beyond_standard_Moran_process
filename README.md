{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# beyond_standard_Moran_process #\n",
    "\n",
    "Here, I am simulating the modified version of well studied, **Moran process (https://en.wikipedia.org/wiki/Moran_process).**\n",
    "\n",
    "**Model:** Consider a haploid population of size $N,$ which is composed of two alleles $A$ and $a$ present in 1 and $N - 1$ copies respectively. The fitnesses of $A$ and $a$ are $1 + s$ and 1 respectively, where $s > 0$. \n",
    "\n",
    "At each discrete time step, one individual is chosen proportional to its fitness from the population.\n",
    "It produces $U$ number of progenies, drawn from the offspring distribution $P_{U}$. In order to keep the population size constant in every time step, $U$ individuals are selected to die randomly from other N − 1 individuals.\n",
    "\n",
    "\n",
    "The offspring distribution ($P_{U}(u)$) is\n",
    "             \\begin{equation}\n",
    "              P_{U}(u) = \\frac{A}{u^{\\gamma}}, \n",
    "              \\end{equation}\n",
    "               where $A$ is the normalization constant.\n",
    "              \n",
    "*In the Moran process, $U$ can take only one value i.e., 2, while in this model, $U$ is a random variable that can take any value ranging from 2 to $N$.*\n",
    "\n",
    "**Important: $\\gamma$ is the parameter of the model. When $\\gamma > 3, $ this model reduces to the standard Moran process.**\n",
    "\n",
    "\n",
    "## Aim: to calculate the probability with which the allele $A$ will take over the population, i.e., the fixation probability of $A$.##\n",
    "\n",
    "![myimage-alt-tag](1.png)\n",
    "\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.6"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
