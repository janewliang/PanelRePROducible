Code to run simulation studies with very large simulated families (families where a large amount of family history is available). Each counselee was simulated to have four maternal aunts, maternal uncles, paternal aunts, paternal uncles, sisters, and brothers. The counselee and their siblings were then simulated to each have four daughters and four sons, for a total of 112 simulated relatives per family. 

- `pp5/`: Simulate families and evaluate them using PanelPRO-5BC, without incorporating risk modifiers.  
- `pp11/`: Simulate families and evaluate them using PanelPRO-11, without incorporating risk modifiers. 

Each sub-directory contains the following files: 
- `fam.R`: Code to simulate 1000 families at a time and run the PanelPRO models. 
- `rscript_fam.job` and`submitJobs_fam.sh`: Shell scripts for running models on the simulated families, for usage on a high performance computing cluster. 
- `diagnostics.sh`: Shell script for obtaining diagnostic metrics, for usage on a high performance computing cluster. 
- `rscript_diagnostics_boot_job`, `submitJobs_diagnostics_boot.sh`, and `combine_diagnostics_boot.sh`: Shell scripts for obtaining diagnostic metrics from 1000 bootstrapped samples, for usage on a high performance computing cluster. 
- `subset_diagnostics.sh`: Shell script for obtaining diagnostic metrics on a subset of the simulated families of the same size as the HCP cohort, for usage on a high performance computing cluster. 
- `rscript_subset_diagnostics_boot_job`, `submitJobs_subset_diagnostics_boot.sh`, and `combine_subset_diagnostics_boot.sh`: Shell scripts for obtaining diagnostic metrics from 1000 bootstrapped samples of the subset of simulated families, for usage on a high performance computing cluster.