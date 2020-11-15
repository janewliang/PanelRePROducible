library(abind)

# Helper functions
source("../../_deps/diagnostics_functions.R")

# Load results from HCP data, to get the counts of each mutation carrier type
load("../../hcp/pp5/results/output/all_dfs.rData")
hcp_prob_df = prob_df 
hcp_mut_df = mut_df

# Load full simulation results
load("results/output/all_dfs.rData")

###############################################################################

# Data frame summarizing counts of each mutation carrier type in the HCP data
hcp_mut_count_df = data.frame(mut=names(hcp_mut_df), 
                              count=colSums(hcp_mut_df))
# Stratified sample of indices based on HCP mutation counts
set.seed(1)
idx = sort(c(
  # Stratified sample of mutation carriers
  unlist(lapply(1:nrow(hcp_mut_count_df), function(i){
    sample(which(mut_df[hcp_mut_count_df$mut[i]]==1), 
           hcp_mut_count_df$count[i], replace=FALSE)
  })), 
  # Sample non-carriers
  sample(which(apply(mut_df==0, 1, all)), 
         nrow(hcp_mut_df)-sum(hcp_mut_count_df$count), replace=FALSE)))

# Subset probabilities and proband mutation information
prob_df = prob_df[,,idx]
mut_df = mut_df[idx,]

save(prob_df, mut_df, file="results/output/subset_all_dfs.rData")

###############################################################################

# Get diagnostics for subsetted data, but for each mutation, drop the 
# carriers for other mutations in the calculations (i.e. only use carriers 
# of mutation of interest and those who are noncarriers for all mutations)

# Indices of noncarriers
noncarrier_idx = rowSums(mut_df)==0

# Put all diagnostics results together
all_out = list(
  ATM = get_all_diagnostics("ATM", NULL, mut_df, prob_df, noncarrier_idx), 
  BRCA1 = get_all_diagnostics("BRCA1", "brcapro", mut_df, prob_df, noncarrier_idx), 
  BRCA2 = get_all_diagnostics("BRCA2", "brcapro", mut_df, prob_df, noncarrier_idx), 
  brcapro_genes = get_all_diagnostics(c("BRCA1", "BRCA2"), "brcapro", 
                                      mut_df, prob_df, noncarrier_idx), 
  CHEK2 = get_all_diagnostics("CHEK2", NULL, mut_df, prob_df, noncarrier_idx), 
  PALB2 = get_all_diagnostics("PALB2", NULL, mut_df, prob_df, noncarrier_idx), 
  Any = get_all_diagnostics("Any", NULL, mut_df, prob_df, noncarrier_idx))

diagnostics = lapply(all_out, function(x){x[[1]]})
carrier_probs = lapply(all_out, function(x){x[2:length(x)]})

# Save output
save(diagnostics, file="results/diagnostics/subset_diagnostics.rData")
save(carrier_probs, file="results/diagnostics/subset_carrier_probs.rData")
