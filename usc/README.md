Code to validate PanelPRO models on data from the USC-Stanford Hereditary Cancer Panel (HCP) Testing Study<sup>[1](#myfootnote1)</sup>. 

### Data processing
- `usc_data/`: Process the HCP data by dropping families where the counselee is a carrier of one or more variants, but not of variants for any of the genes in the model; or where the counselee has a variant of uncertain significance (VUS), but not a pathogenic variant for any of the genes in the model. Format the data so that it is ready to be passed into the models. 
  - `pp5/`: Contains the file `get_usc_data.R`, which is used to process the HCP data for validation with PanelPRO-5BC. 
  - `pp7/`: Contains the file `get_usc_data.R`, which is used to process the HCP data for validation with PanelPRO-7. 
  - `pp11/`: Contains the file `get_usc_data.R`, which is used to process the HCP data for validation with PanelPRO-11. 
  - `usc_dat.RData`: Data for the HCP cohort of 2000 families, which includes detailed pedigrees, family history of cancer, and panel gene testing results. Not publicly available. 
  - `usc_summary.R`: Code to obtain summary information about the HCP cohort, such as the number of probands who are carriers of pathogenic variants, the number of probands who are cancer cases, and the number of relatives with tumor markers. 

### Model validation
- `pp5/`: Validate PanelPRO-5BC on the HCP cohort, without incorporating risk modifiers. 
- `ppt_riskmod/`: Validate PanelPRO-5BC on the HCP cohort, while incorporating tumor marker testing as an additional risk modifier. 
- `pp11/`: Validate PanelPRO-11 on the HCP cohort, without incorporating risk modifiers. 
- `pp11_riskmod/`: Validate PanelPRO-11 on the HCP cohort, while incorporating tumor marker testing as an additional risk modifier. 

Each sub-directory contains the following files: 
- `fam.R`: Code to split the HCP pedigrees into 20 approximately equal-sized parts and run the PanelPRO models. 
- `rscript_fam.job` and`submitJobs_fam.sh`: Shell scripts for running models on the HCP cohort, for usage on a high performance computing cluster. 
- `diagnostics.sh`: Shell script for obtaining diagnostic metrics, for usage on a high performance computing cluster. 
- `rscript_diagnostics_boot_job`, `submitJobs_diagnostics_boot.sh`, and `combine_diagnostics_boot.sh`: Shell scripts for obtaining diagnostic metrics from 1000 bootstrap samples, for usage on a high performance computing cluster. 

---

<a name="myfootnote1">1</a>. Idos, G. E., Kurian, A. W., Ricker, C., Sturgeon, D., Culver, J. O., Kingham, K. E., ... & Levonian, P. (2019). Multicenter prospective cohort study of the diagnostic yield and patient experience of multiplex gene panel testing for hereditary cancer risk. JCO Precision Oncology, 3, 1-12.