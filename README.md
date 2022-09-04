# est-se

This repository contains the data and code used for the analyses in “Standard error estimation in meta-analysis of studies reporting medians”. 

The analyses were originally run on R version 4.1.3 with packages ``data.table`` (version 1.14.2), ``doParallel`` (version 1.0.17), ``estmeansd`` (version 0.2.1), ``fdrtool`` (version 1.2.17), ``foreach`` (version 1.5.2), ``metaBLUE`` (version 1.0.0), ``metafor`` (version 1.6-0), ``stringr`` (version 1.7.6), and ``xtable`` (version 1.8-4). 

The study-level simulations took approximately 6 hours to run when parallelized across 20 CPU cores. The meta-analysis simulations took approximately 26 hours to run when parallelized across 40 CPU cores. The analyses for obtaining the true $I^2$ values took approximate 3 hours to run when parallelized across 10 CPU cores. The data application took approximately 15 minutes to run when parallelized across 9 CPU cores.

## Simulation Study

Before running the simulations, one must first run the ``settings.R`` file which saves the simulation settings in the ``simulation_settings.RData`` file. The study-level simulations can be run by running the ``main_studylevel.R`` file. The meta-analysis simulations can be run by running the ``main_ma.R`` file. The true values of $I^2$ can be obtained by running the ``true_I2.R`` file. All simulation results can be analyzed by running the ``analyze_results_simulations.R`` file.

The following helper files are used:

-   ``helper_functions_simulation.R``: Contains helper functions used in the simulation study (e.g., applying the transformation-based approaches, simulating data, and analyzing simulation results)
-   ``mln.R``: Implements the MLN transformation-based approach. Thanks to Siyu Cai for providing much of this code.

## Data Application

The data cleaning and application of the transformation-based methods is performed in the ``main_application.R`` file. This file uses helper functions in the ``helper_functions_application.R``file. One can generate the figures and tables corresponding to the data application by running the ``analyze_results_application.R`` file. 
