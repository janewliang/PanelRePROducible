library(abind)
library(jsmodule)

# Helper functions
source("../../_deps/diagnostics_functions.R")

###############################################################################

# 1000 sets of family simulations
K = 1000

# Putting all of the family simulation results together
all_fam_output = NULL
for (k in 1:K){
  load(paste0("results/output/fam", k, ".rData"))
  all_fam_output = c(all_fam_output, fam_output)
}
# Simulated families
usc_families_sim = lapply(all_fam_output, function(x){ x$fam})
# Carrier probabilities
all_probs = lapply(all_fam_output, function(x){ x$probs})

# No failed families
failed_families = sapply(all_probs, function(x){length(x)==0 || class(x)=="try-error"})
all_probs[failed_families] = NULL

# Pull out the proband mutation information from each family
mutations = c("ATM", "BRCA1", "BRCA2", "CHEK2", "PALB2")
mut_df = data.frame(t(sapply(usc_families_sim, function(fam){
  unlist(fam[fam$isProband==1,mutations])
})))
mut_df = mut_df[!failed_families,]

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
  brcapro_genes = get_all_diagnostics(c("BRCA1", "BRCA2"), "brcapro", 
                                      mut_df, prob_df, noncarrier_idx), 
  CHEK2 = get_all_diagnostics("CHEK2", NULL, mut_df, prob_df, noncarrier_idx), 
  PALB2 = get_all_diagnostics("PALB2", NULL, mut_df, prob_df, noncarrier_idx), 
  Any = get_all_diagnostics("Any", NULL, mut_df, prob_df, noncarrier_idx))

diagnostics = lapply(all_out, function(x){x[[1]]})
carrier_probs = lapply(all_out, function(x){x[2:length(x)]})

# BayesMendel genes
BMgenes = c("BRCA1", "BRCA2", "brcapro_genes")
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
save(diagnostics, netreclass, file="results/diagnostics/diagnostics.rData")
save(carrier_probs, file="results/diagnostics/carrier_probs.rData")
