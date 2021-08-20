library(abind)

# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Helper functions
source("../../_deps/diagnostics_functions.R")

# Load results
load("results/output/all_dfs.rData")

###############################################################################

# Bootstrap diagnostics
R = 10 # Number of bootstrap replicates
n = nrow(mut_df) # Number of probands/families

# Get diagnostics from bootstrapped samples
diagnostics_boot = lapply(1:R, function(i){
  # Set random seed
  set.seed(s*R+i)
  
  # Sample observations with replacement
  idx = sample(n, n, replace=TRUE)
  
  # Bootstrap data
  boot_mut_df = mut_df[idx,]
  boot_prob_df = prob_df[,,idx]
  
  # Indices of noncarriers
  noncarrier_idx = rowSums(boot_mut_df)==0
  
  # Get bootstrap diagnostics
  all_out = abind(
    ATM = get_all_diagnostics("ATM", NULL, boot_mut_df, boot_prob_df, 
                              noncarrier_idx, return_carriers=FALSE), 
    BRCA1 = get_all_diagnostics("BRCA1", "brcapro", boot_mut_df, boot_prob_df, 
                                noncarrier_idx, return_carriers=FALSE), 
    BRCA2 = get_all_diagnostics("BRCA2", "brcapro", boot_mut_df, boot_prob_df, 
                                noncarrier_idx, return_carriers=FALSE), 
    BRCAPRO_genes = get_all_diagnostics(c("BRCA1", "BRCA2"), "brcapro", 
                                        boot_mut_df, boot_prob_df, 
                                        noncarrier_idx, return_carriers=FALSE), 
    CDKN2A = get_all_diagnostics("CDKN2A", "melapro", boot_mut_df, boot_prob_df, 
                                 noncarrier_idx, return_carriers=FALSE), 
    CHEK2 = get_all_diagnostics("CHEK2", NULL, boot_mut_df, boot_prob_df, 
                                noncarrier_idx, return_carriers=FALSE), 
    EPCAM = get_all_diagnostics("EPCAM", NULL, boot_mut_df, boot_prob_df, 
                                noncarrier_idx, return_carriers=FALSE), 
    MLH1 = get_all_diagnostics("MLH1", "mmrpro", boot_mut_df, boot_prob_df, 
                               noncarrier_idx, return_carriers=FALSE),
    MSH2 = get_all_diagnostics("MSH2", "mmrpro", boot_mut_df, boot_prob_df, 
                               noncarrier_idx, return_carriers=FALSE),
    MSH6 = get_all_diagnostics("MSH6","mmrpro", boot_mut_df, boot_prob_df, 
                               noncarrier_idx, return_carriers=FALSE), 
    MMRpro_genes = get_all_diagnostics(c("MLH1", "MSH2", "MSH6"), "mmrpro", 
                                       boot_mut_df, boot_prob_df, 
                                       noncarrier_idx, return_carriers=FALSE),
    PALB2 = get_all_diagnostics("PALB2", NULL, boot_mut_df, boot_prob_df, 
                                noncarrier_idx, return_carriers=FALSE), 
    PMS2 = get_all_diagnostics("PMS2", NULL, boot_mut_df, boot_prob_df, 
                               noncarrier_idx, return_carriers=FALSE), 
    Any_BM = get_all_diagnostics("Any_BM", c("brcapro", "mmrpro", "melapro"), 
                                 boot_mut_df, boot_prob_df, 
                                 noncarrier_idx, return_carriers=FALSE), 
    Any = get_all_diagnostics("Any", NULL, boot_mut_df, boot_prob_df, 
                              noncarrier_idx, return_carriers=FALSE), 
    along=3) 
  return(all_out)
})

save(diagnostics_boot, file = paste0("results/diagnostics/diagnostics_boot", 
                                     s, ".rData"))
