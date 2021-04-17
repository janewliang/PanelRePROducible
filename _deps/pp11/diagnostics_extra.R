library(abind)

# Helper functions
source("../../_deps/diagnostics_functions.R")

###############################################################################

# 1000 sets of family simulations
K = 1000

# Putting all of the family simulation results together
all_fam_output = NULL
# Original 1000
for (k in 1:K){
  load(paste0("../pp11/results/output/fam", k, ".rData"))
  all_fam_output = c(all_fam_output, fam_output)
}
# Extra 1000
for (k in 1:K){
  load(paste0("results/output/fam", k, ".rData"))
  all_fam_output = c(all_fam_output, fam_output)
}
# Simulated families
usc_families_sim = lapply(all_fam_output, function(x){ x$fam})
# Carrier probabilities
all_probs = lapply(all_fam_output, function(x){ x$probs})

# Pull out the proband mutation information from each family
mutations = c("ATM", "BRCA1", "BRCA2", "CDKN2A", "CHEK2", "EPCAM",
              "MLH1", "MSH2", "MSH6", "PALB2", "PMS2")
mut_df = data.frame(t(sapply(usc_families_sim, function(fam){
  unlist(fam[fam$isProband==1,mutations])
})))

# Pull out the probabilities for each model and save as a 3D array
prob_df = abind(all_probs, along=3)

save(prob_df, mut_df, file="results/output/all_dfs.rData")

###############################################################################

# Get diagnostics, but for each mutation, drop the carriers for other mutations
# in the calculations (i.e. only use carriers of mutation of interest and those 
# who are noncarriers for all mutations)

# Indices of noncarriers
noncarrier_idx = rowSums(mut_df)==0

# Put all diagnostics results together
all_out = list(
  ATM = get_all_diagnostics("ATM", NULL, mut_df, prob_df, noncarrier_idx), 
  BRCA1 = get_all_diagnostics("BRCA1", "brcapro", mut_df, prob_df, noncarrier_idx), 
  BRCA2 = get_all_diagnostics("BRCA2", "brcapro", mut_df, prob_df, noncarrier_idx), 
  BRCAPRO_genes = get_all_diagnostics(c("BRCA1", "BRCA2"), "brcapro", 
                                      mut_df, prob_df, noncarrier_idx), 
  CDKN2A = get_all_diagnostics("CDKN2A", "melapro", mut_df, prob_df, noncarrier_idx), 
  CHEK2 = get_all_diagnostics("CHEK2", NULL, mut_df, prob_df, noncarrier_idx), 
  EPCAM = get_all_diagnostics("EPCAM", NULL, mut_df, prob_df, noncarrier_idx), 
  MLH1 = get_all_diagnostics("MLH1", "mmrpro", mut_df, prob_df, noncarrier_idx),
  MSH2 = get_all_diagnostics("MSH2", "mmrpro", mut_df, prob_df, noncarrier_idx),
  MSH6 = get_all_diagnostics("MSH6","mmrpro", mut_df, prob_df, noncarrier_idx), 
  MMRpro_genes = get_all_diagnostics(c("MLH1", "MSH2", "MSH6"), "mmrpro", 
                                     mut_df, prob_df, noncarrier_idx),
  PALB2 = get_all_diagnostics("PALB2", NULL, mut_df, prob_df, noncarrier_idx), 
  PMS2 = get_all_diagnostics("PMS2", NULL, mut_df, prob_df, noncarrier_idx), 
  Any_BM = get_all_diagnostics("Any_BM", c("brcapro", "mmrpro", "melapro"), 
                               mut_df, prob_df, noncarrier_idx), 
  Any = get_all_diagnostics("Any", NULL, mut_df, prob_df, noncarrier_idx)) 

diagnostics = lapply(all_out, function(x){x[[1]]})
carrier_probs = lapply(all_out, function(x){x[2:length(x)]})

# Save output
save(diagnostics, file="results/diagnostics/diagnostics.rData")
save(carrier_probs, file="results/diagnostics/carrier_probs.rData")

###############################################################################

# Get carrier counts
mut_counts = colSums(mut_df)

# Get carrier counts for each bootstrap
n = nrow(mut_df) # Number of probands/families
boot_mut_counts = t(sapply(1:1000, function(i){
  
  # Set random seed
  set.seed(i+n)
  
  # Sample observations with replacement
  idx = sample(n, n, replace=TRUE)
  
  # Bootstrap data
  return(colSums(mut_df[idx,]))
}))
save(mut_counts, boot_mut_counts, 
     file = "results/diagnostics/boot_mut_counts.rData")
