# KTLC (Koopmans’ Theorem-Compliant Long-range Corrected) density functional

KTLC is a new density functional that uses the optimal range-separated parameter μ for a given molecule [1]. The optimal μ can be easily predicted using [Google colaboratory (ColabKTLC)](https://colab.research.google.com/github/molecule-generator-collection/KTLC/blob/main/ColabKTLC.ipynb). 

[1] K. Terayama, Y. Osaki, T. Fujita, R. Tamura, M. Naito, K. Tsuda, T. Matsui, M. Sumita, Koopmans’ Theorem-Compliant Long-range Corrected (KTLC) Density Functional Mediated by Black-box Optimization and Machine Learning for Organic Molecules, arXiv, 2022

# Optimization of the range-separated parameters using Bayesian opitmization

The optimisation of the range-separated parameters for a given molecule using Bayesian optimisation can be performed using param_opt_main.py as follows.

## Requirements

1. [Gaussian](https://gaussian.com)==16
2. [Python](https://www.anaconda.com/download/)==3.7
3. [rdkit-pypi](https://anaconda.org/rdkit/rdkit)==2021.03.5
4. [QCforever](https://github.com/molecule-generator-collection/QCforever)
5. [chem_LAQA](https://github.com/inter-info-lab/chem_laqa)
6. [xtb](https://github.com/grimme-lab/xtb)==6.4.1
7. [bayesian-optimization](https://github.com/bayesian-optimization/BayesianOptimization)==1.4.3
