# PanelRePROducible

This repository contains the code to reproduce the results presented in the paper 

  Liang, J. W., Idos, G., Hong, C., Gruber, S., Parmigiani, G., & Braun, D. (2020). PanelPRO: a general framework for multi-gene, multi-cancer Mendelian risk prediction models. 

## Dependencies
- `PanelPRO` R package, Lee, G., et al. (2020)<sup>[1](#myfootnote1)</sup>.
- `BayesMendel` R package, Chen, S., et al. (2004)<sup>[2](#myfootnote2)</sup>. Available for download [here](https://projects.iq.harvard.edu/bayesmendel/bayesmendel-r-package).
- Other R packages: [`abind`](https://cran.r-project.org/web/packages/abind/index.html), [`knitr`](https://cran.r-project.org/web/packages/knitr/index.html), [`pROC`](https://cran.r-project.org/web/packages/pROC/index.html)
- Data from the USC-Stanford Hereditary Cancer Panel (HCP) Testing Study. Idos, G., et al. (2019)<sup>[3](#myfootnote3)</sup>. Not publicly available.

## Navigation
Additional details can be found in the sub-directory README files. 

### Shared Code Dependencies
- `_deps/`: Functions for simulating pedigrees, running PanelPRO models, and obtaining diagnostic metrics. They are called as dependencies to run the simulation studies and data application. 

### Simulation Studies
- `sim/`: Code to run simulation studies where the simulated families have family structures sampled from the HCP cohort. 
- `sim_upper/`: Code to run simulation studies with very large simulated families (families where a large amount of family history is available). 

### Data Application
- `hcp/`: Code to validate PanelPRO models on the HCP cohort. 

### Visualizations
- `plots/`: Code to generate figures and tables in paper based on results from `sim/`, `sim_upper`, and `hcp`. 

## Useful Terminology
- *PanelPRO-5BC*: Mendelian risk prediction model that identifies individuals at high risk for breast and ovarian cancer due to ATM, BRCA1, BRCA2, CHEK2, and PALB2 mutations. 
- *PanelPRO-7*: Mendelian risk prediction model that identifies individuals at high risk for breast, colorectal, endometrial, ovarian, pancreatic, and skin cancer due to BRCA1, BRCA2, CDKN2A, MLH1, MSH2, MSH6, and (hypothetical) PANC mutations. This model spans all of the genes and cancers currently incorporated in models in the BayesMendel R package<sup>[2](#myfootnote2)</sup>. 
- *PanelPRO-11*: Mendelian risk prediction model that identifies individuals at high risk for eleven cancers (brain, breast, colorectal, endometrial, gastric, kidney, melanoma, ovarian, pancreatic, prostate, and small intestine) due to mutations in eleven genes (ATM, BRCA1, BRCA2, CDKN2A, CHEK2, EPCAM, MSH1, MSH2, MLH6, PALB2, and PMS2). 
- *Diagnostic metrics*: To evaluate the posterior carrier probabilities from our models, we used the area under the curve (AUC) as a measure for discrimination, the expected divided by the observed number of events as a measure of calibration, and mean squared error (MSE) as a measure of accuracy. When computing these diagnostics, non-carriers were defined as individuals who are not carriers of any mutation in the model. 

---

<a name="myfootnote1">1</a>. Lee, G., Zhang, Q., Liang, J. W., Huang, T., Choirat, C., Parmigiani, G., & Braun, D. (2020). PanelPRO: A R package for multi-syndrome, multi-gene risk modeling for individuals with a family history of cancer. arXiv preprint arXiv:2010.13011.

<a name="myfootnote2">2</a>. Chen, S., Wang, W., Broman, K. W., Katki, H. A., & Parmigiani, G. (2004). BayesMendel: an R environment for Mendelian risk prediction. Statistical applications in genetics and molecular biology, 3(1).

<a name="myfootnote3">3</a>. Idos, G. E., Kurian, A. W., Ricker, C., Sturgeon, D., Culver, J. O., Kingham, K. E., ... & Levonian, P. (2019). Multicenter prospective cohort study of the diagnostic yield and patient experience of multiplex gene panel testing for hereditary cancer risk. JCO Precision Oncology, 3, 1-12.