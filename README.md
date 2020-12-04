# tumorgrowth - Analysis of the equilibrium phase in immune-controlled tumors to predict best strategies for cancer treatment

Source codes to compute and perform sensitivity analysis of the stationary profile in a size and space structured tumor-immune interaction model.

# Program requirements
the following python version must be installed to run to source codes:
- python 3.7

# Python packages 
- fipy (version 3.3)
- gmsh (version 3.0.6)
- scipy (version 1.5.3)
- matplotlib (version 3.2.1)
- numpy (version 1.19.1)
- pandas (version 1.1.0)
- seaborn (version 0.10.1)
- h5py (version 2.10.0)
- pygpc (version 0.2.7.5)
- cma (version 3.0.3)

# The organization of the code

1. Simulations of the evolution problem 
- the file `2DTumorImmuneInteraction.py` 
- the file `2DTumorImmuneInteraction_withRealParameters.py`

2. Estimation the model parameters
- the file `params_estimation.py`

3. The power-dichotomy algorithm
- the file `power_dichotomy_algorithm.py`
- the file `power_dichotomy_algorithm_withRealParameters.py`
- the file `power_dichotomy_algorithm_withRealParameters_multSimu.py`
- the file `power_dichotomy_algorithm_withRealParameters_multSimu_matrix.py`

4. The global sensitivity analysis of the tumor mass at equilibrium
- the file `sensitivity_analysis_withRealParameters.py`
- the file `sensitivity_analysis_results_postprocessing.py`
