# KTLC (Koopmans’ Theorem-Compliant Long-range Corrected) density functional

KTLC is a new density functional that uses the optimal range-separated parameter μ for a given molecule [1]. The optimal μ can be easily predicted using [Google colaboratory (ColabKTLC)](https://colab.research.google.com/github/molecule-generator-collection/KTLC/blob/main/ColabKTLC.ipynb). 

[1] K. Terayama, Y. Osaki, T. Fujita, R. Tamura, M. Naito, K. Tsuda, T. Matsui, M. Sumita, [Koopmans’ Theorem-Compliant Long-range Corrected (KTLC) Density Functional Mediated by Black-box Optimization and Data-Driven Prediction for Organic Molecules](https://doi.org/10.26434/chemrxiv-2023-5nktj), ChemRxiv, 2022. [DOI: 10.26434/chemrxiv-2023-5nktj]

# Optimization of the range-separated parameters using Bayesian opitmization

The optimisation of the range-separated parameters for a given molecule using Bayesian optimisation can be performed using param_opt_main.py as follows.

## Requirements

1. [Gaussian](https://gaussian.com)==16
2. [Python](https://www.anaconda.com/download/)==3.7
3. [rdkit-pypi](https://anaconda.org/rdkit/rdkit)==2021.03.5
4. [QCforever](https://github.com/molecule-generator-collection/QCforever)==1.0.0
5. [chem LAQA](https://github.com/inter-info-lab/chem_laqa)
6. [xtb](https://github.com/grimme-lab/xtb)==6.4.1
7. [bayesian-optimization](https://github.com/bayesian-optimization/BayesianOptimization)==1.4.3
8. [openbabel](https://github.com/openbabel/openbabel) with eigen
9. [pyyaml](https://github.com/yaml/pyyaml) == 6.0

##  Preparation of config file

See an example of the config file `config/setting_test.yaml`. 

## Optimization of the range-separated parameters

```bash
python param_opt_main -c config/setting_test.yaml
```

## Result of the optimization

All results are stored in the \{result_prefix\}\_result\{index\} directory. Here, "result_prefix" and "index" are the options in the config file. 
The summaries of the structural optimization, the optimization of the range-separated parameters, and the final result with the optimized parameters are contained in init_calc_result_dict.pkl, bo_calc_result_dict.pkl, and best_param_calc_result_dict.pkl, respectively. 


## Support option in the config file

|Option|Description|Option in the paper| 
|---|---|---|
|`index`|Index for a target molecule|-|
|`smiles`|SMILES of a target molecule|-|
|`init_Gau_option`|Options of the conformational optimization using QCforever. See the details of QCforever.|opt satkoopmans deen vip vea homolumo dipole uv|
|`init_functional`|Functional for the conformational optimization. Typical options are B3LYP and BLYP. See QCforever for the details.|B3LYP|
|`init_basis`|Basis set for the conformational optimization. Typical options are 6-31G*, 3-21G*, and STO3G.|6-31G*|
|`init_solvent`|Solvent for the conformational optimization. "0" indicates water.|0|
|`opt_Gau_option`|Options of the optimization of the the range-separated parameters. "satkoopmans" is requred for this optimization.|satkoopmans homolumo|
|`opt_functional`|Functional for the paramter optimization. Typecal options are LC-BLYP and CAM-B3LYP|LC-BLYP or CAM-B3LYP|
|`opt_basis`|Basis set for the the paramter optimization. |6-31G*|
|`opt_solvent`|Solvent for the paramter optimization.|0|
|`calc_best_parameter`|If True, the calcaulation with the optimized parameters are performed using the "best_Gau_option".|True|
|`best_Gau_option`|Options for the calculation with the optimized parameters.|satkoopmans deen vip vea homolumo dipole uv|
|`bo_strategy`|Strategy of Bayesian optimization. EI and random are implemented.|EI or random|
|`bo_init_num`|Number of initialization for Bayesian optimization.|3 (LC) or 4 (CAM)|
|`bo_search_num`|Number of iterations in Bayesian optimization.|50 (LC) or 100 (CAM)|
|`bo_stop_threshold`|Threshold (eV) for stopping the optiomization. |0.01|

