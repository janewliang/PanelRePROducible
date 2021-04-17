Functions for running PanelPRO models and obtaining diagnostic metrics. They are called as dependencies to run the simulation studies and data application. 

### Functions and Code Shared Across Models
- `BackCompatibleDatabase.R`: Code to generate databases of model parameters/inputs that will return from PanelPRO<sup>[1](#myfootnote1)</sup> models that are back-compatible with models from the BayesMendel package<sup>[2](#myfootnote2)</sup>. 
- `run_model_functions.R`:  Functions that run models and extract carrier probability results for the simulation studies and data application. 
- `diagnostics_functions.R`: Functions for obtaining diagnostic metrics. 
- `combine_diagnostics_boot.R`: Code to load bootstrapped diagnostic metrics returned from separate cluster jobs. 
- `combine_diagnostics_boot_extra.R`: Code to load bootstrapped diagnostic metrics returned from separate cluster jobs, specifically for the results from running extra PanelPRO-11 simulations. 
- `combine_subset_diagnostics_boot.R`: Code to load bootstrapped diagnostic metrics for a subset of families returned from separate cluster jobs. 

### Model-Specific Functions and Code
- `pp5/`: Functions and code specific to running PanelPRO-5BC. 
- `pp7/`: Functions and code specific to running PanelPRO-7. 
- `pp11/`: Functions and code specific to running PanelPRO-11. 

Each sub-directory contains the following files: 
- `diagnostics.R`: Code to obtain diagnostic metrics from evaluating models on simulated families. 
- `diagnostics_boot.R`: Code to obtain diagnostic metrics from 10 bootstrap samples. 

In addition, `pp5/` and `pp11/` contain the following files: 
- `diagnostics_usc.R`: Code to obtain diagnostic metrics from evaluating models on families from the HCP cohort. 
- `subset_diagnostics.R`: Code to obtain diagnostic metrics from evaluating models on a subset of simulated families. 
- `subset_diagnostics_boot.R`: Code to obtain diagnostic metrics from 10 bootstrap samples of a subset of simulated families. 

In addition, `pp11/` contains the following files: 
- `diagnostics_extra.R`: Code to obtain diagnostic metrics from evaluating models on simulated families, for all 2 million families simulated and evaluated under PanelPRO-11. 
- `diagnostics_boot_extra.R`: Code to obtain diagnostic metrics from 10 bootstrap samples. Bootstrap samples are taken from an original sample size ranging from 1 million to 2 million in 200,000 increments. 

---

<a name="myfootnote1">1</a>. Lee, G., Zhang, Q., Liang, J. W., Huang, T., Choirat, C., Parmigiani, G., & Braun, D. (2020). PanelPRO: A R package for multi-syndrome, multi-gene risk modeling for individuals with a family history of cancer. arXiv preprint arXiv:2010.13011.

<a name="myfootnote2">2</a>. Chen, S., Wang, W., Broman, K. W., Katki, H. A., & Parmigiani, G. (2004). BayesMendel: an R environment for Mendelian risk prediction. Statistical applications in genetics and molecular biology, 3(1).