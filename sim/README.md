Code to run simulation studies where the simulated families have family structures sampled from the HCP cohort. 

- `pp5/`: Simulate families and evaluate them using PanelPRO-5BC, without incorporating risk modifiers. 
- `ppt_riskmod/`: Simulate families and evaluate them using PanelPRO-5BC, while incorporating tumor marker testing as an additional risk modifier. 
- `pp7_gen_full`: Simulate families (under assumptions of the full PanelPRO-7 model) and evaluate them using PanelPRO-7, without incorporating risk modifiers. 
- `pp7_gen_sub`: Simulate families (under assumptions of the sub-models nested in PanelPRO-7) and evaluate them using PanelPRO-7, without incorporating risk modifiers. 
- `pp11/`: Simulate families and evaluate them using PanelPRO-11, without incorporating risk modifiers. 
- `pp11_riskmod/`: Simulate families and evaluate them using PanelPRO-11, while incorporating tumor marker testing as an additional risk modifier. 

Each sub-directory contains the following files: 
- `fam.R`: Code to simulate 1000 families at a time and run the PanelPRO models. 
- `rscript_fam.job` and`submitJobs_fam.sh`: Shell scripts for running models on the simulated families, for usage on a high performance computing cluster. 
- `diagnostics.sh`: Shell script for obtaining diagnostic metrics, for usage on a high performance computing cluster. 
- `rscript_diagnostics_boot_job`, `submitJobs_diagnostics_boot.sh`, and `combine_diagnostics_boot.sh`: Shell scripts for obtaining diagnostic metrics from 1000 bootstrapped samples, for usage on a high performance computing cluster. 

In addition, `pp5/` and `pp11/` contain the following files: 
- `subset_diagnostics.sh`: Shell script for obtaining diagnostic metrics on a subset of the simulated families of the same size as the HCP cohort, for usage on a high performance computing cluster. 
- `rscript_subset_diagnostics_boot_job`, `submitJobs_subset_diagnostics_boot.sh`, and `combine_subset_diagnostics_boot.sh`: Shell scripts for obtaining diagnostic metrics from 1000 bootstrapped samples of the subset of simulated families, for usage on a high performance computing cluster. 