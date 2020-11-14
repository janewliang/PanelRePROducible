# PanelRePROducible

This repository contains the code to reproduce the results presented in the paper "PanelPRO: a general framework for multi-gene, multi-cancer Mendelian risk prediction models". 

## Dependencies
- `PanelPRO` R package [^fn1]
- `BayesMendel` R package [^fn2], available [here](https://projects.iq.harvard.edu/bayesmendel/bayesmendel-r-package) 
- Other R packages: `abind`, `jsmodule`, `knitr`, `pROC`
- Data from the USC-Stanford Hereditary Cancer Panel (HCP) Testing Study [^fn3], not publicly  available

## Navigation
Additional details can be found in the sub-directory READMEs. 

### Shared Functions
- `_deps/`: Functions for running models and obtaining model diagnostics, shared across different models, simulations, and data applications. 
- `_scripts/`: Functions for obtaining and summarizing model diagnostics, specific to each type of model, simulation study, and data application. 

### Simulations
- `simulate_families`: Code to simulate pedigrees, including genotypes, family history of cancer, and tumor marker testing results. 
- `sim/`: Code to run simulation studies with family structures sampled from the HCP cohort. 
- `big_fam/`: Code to run simulation studies with very large families. 

### Data Application
- `usc/`: Code to validate models on the HCP cohort. 

### Visualizations
- `plots/`: Code to generate figures and tables in paper. 

---

[^fn1]: Lee, G., Zhang, Q., Liang, J. W., Huang, T., Choirat, C., Parmigiani, G., & Braun, D. (2020). PanelPRO: A R package for multi-syndrome, multi-gene risk modeling for individuals with a family history of cancer. arXiv preprint arXiv:2010.13011.

[^fn2]: Chen, S., Wang, W., Broman, K. W., Katki, H. A., & Parmigiani, G. (2004). BayesMendel: an R environment for Mendelian risk prediction. Statistical applications in genetics and molecular biology, 3(1).

[^fn3]: Idos, G. E., Kurian, A. W., Ricker, C., Sturgeon, D., Culver, J. O., Kingham, K. E., ... & Levonian, P. (2019). Multicenter prospective cohort study of the diagnostic yield and patient experience of multiplex gene panel testing for hereditary cancer risk. JCO Precision Oncology, 3, 1-12.
