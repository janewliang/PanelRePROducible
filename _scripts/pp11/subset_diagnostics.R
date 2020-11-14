library(abind)
library(jsmodule)

# Helper functions
source("../../_deps/diagnostics_functions.R")

# Load results from USC data, to get the counts of each mutation carrier type
load("../../usc/pp11/results/output/all_dfs.rData")
usc_prob_df = prob_df 
usc_mut_df = mut_df

# Load full simulation results
load("results/output/all_dfs.rData")

###############################################################################

# Data frame summarizing counts of each mutation carrier type in the USC data
usc_mut_count_df = data.frame(mut=names(usc_mut_df), 
                              count=colSums(usc_mut_df))
# Stratified sample of indices based on USC mutation counts
set.seed(1)
idx = sort(c(
  # Stratified sample of mutation carriers
  unlist(lapply(1:nrow(usc_mut_count_df), function(i){
    sample(which(mut_df[usc_mut_count_df$mut[i]]==1), 
           usc_mut_count_df$count[i], replace=FALSE)
  })), 
  # Sample non-carriers
  sample(which(apply(mut_df==0, 1, all)), 
         nrow(usc_mut_df)-sum(usc_mut_count_df$count), replace=FALSE)))

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

# BayesMendel genes
BMgenes = c("BRCA1", "BRCA2", "brcapro_genes", 
            "MLH1", "MSH2", "MSH6", "MMRpro_genes", 
            "CDKN2A", "Any_BM")
# Thresholds for being a carrier
thresholds = c("thresh=0.01"=0.01, "thresh=0.025"=0.025)
# Net reclassification
netreclass = lapply(thresholds, function(thresh){
  lapply(BMgenes, function(gene){
    df = data.frame(true=carrier_probs[[gene]]$mutation_status, 
                    pp_sub=carrier_probs[[gene]]$carrier_PP_sub, 
                    pp=carrier_probs[[gene]]$carrier_PP)
    try(reclassificationJS(data=df, cOutcome=1, 
                           predrisk1=df$pp_sub, predrisk2=df$pp, 
                           cutoff=c(0, thresh, 1), dec.value=7, dec.p=7), 
        silent=TRUE)
  })
})

# Save output
save(diagnostics, netreclass, file="results/diagnostics/subset_diagnostics.rData")
save(carrier_probs, file="results/diagnostics/subset_carrier_probs.rData")
