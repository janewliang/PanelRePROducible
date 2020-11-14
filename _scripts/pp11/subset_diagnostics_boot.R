library(abind)

# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Helper functions
source("../../_deps/diagnostics_functions.R")

# Load subsetted simulation results
load("results/output/subset_all_dfs.rData")

###############################################################################

# Bootstrap diagnostics
R = 10 # Number of bootstrap replicates
n = nrow(mut_df) # Number of probands/families

# Get diagnostics from bootstrapped samples of subsetted data
diagnostics_boot = lapply(1:R, function(i){
  # Set random seed
  set.seed(s*R+i)
  
  # Pull out probands that have mutations
  idx1 = apply(mut_df, 2, function(x){
    which(x>0)
  })
  # Find the indices of the probands that have no mutations
  idx0 = which(rowSums(mut_df) == 0)
  
  # Number of cases
  n1 = n - length(idx0)
  # Number of non-cases
  n0 = length(idx0)
  
  # Sample observations with replacement
  # Some people are carriers of multiple mutations, so the bootstrapped sample
  # may be a slightly different size and/or have a very slightly different 
  # number of mutation carriers
  idx = c(unlist(sapply(idx1, function(x){
    sample(x, length(x), replace=TRUE)
  })), sample(idx0, n0, replace=TRUE))
  
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
    brcapro_genes = get_all_diagnostics(c("BRCA1", "BRCA2"), "brcapro", 
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

save(diagnostics_boot, 
     file = paste0("results/diagnostics/subset_diagnostics_boot", 
                   s, ".rData"))
