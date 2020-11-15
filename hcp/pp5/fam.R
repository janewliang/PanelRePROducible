# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Load functions for running models
source("../../_deps/run_model_functions.R")

# Load HCP families
load("../hcp_data/pp5/hcp_families.rData")

###############################################################################

# Split the families for parallelization
idx = split(1:length(hcp_families_PP), 
            ceiling(seq_along(1:length(hcp_families_PP)) / 
                      ceiling(length(hcp_families_PP)/20)))


# Run models for comparison on families
fam_output = lapply(idx[[s]], function(i){
  # Extract family
  fam_PP = hcp_families_PP[[i]]
  fam_BM = hcp_families_BM[[i]]
  
  # Run models
  out = try(run_models_5BC(fam_PP, fam_BM), TRUE)
  return(out)
})


# Save results
save(fam_output, file = paste0("results/output/fam", s, ".rData"))
