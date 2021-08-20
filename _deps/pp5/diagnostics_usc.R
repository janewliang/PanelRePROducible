library(abind)

# Helper functions
source("../../_deps/diagnostics_functions.R")

# Load USC families
load("../usc_data/pp5/usc_families.rData")

###############################################################################

# 233 sets of family results
S = 233

# Putting all of the family simulation results together
all_probs = NULL
for (s in 1:S){
  load(paste0("results/output/fam", s, ".rData"))
  all_probs = c(all_probs, fam_output)
}
# 17 failed families out of 1629
failed_families = sapply(all_probs, function(x){length(x)==0 || class(x)=="try-error"})
failed_errs = all_probs[failed_families]
names(failed_errs) = which(failed_families)
all_probs[failed_families] = NULL

# Pull out the proband mutation information from each family
mut_df = proband_muts[!failed_families,]

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
  CHEK2 = get_all_diagnostics("CHEK2", NULL, mut_df, prob_df, noncarrier_idx), 
  PALB2 = get_all_diagnostics("PALB2", NULL, mut_df, prob_df, noncarrier_idx), 
  Any = get_all_diagnostics("Any", NULL, mut_df, prob_df, noncarrier_idx))

diagnostics = lapply(all_out, function(x){x[[1]]})
carrier_probs = lapply(all_out, function(x){x[2:length(x)]})

# Save output
save(diagnostics, failed_errs, failed_families, 
     file="results/diagnostics/diagnostics.rData")
save(carrier_probs, file="results/diagnostics/carrier_probs.rData")
