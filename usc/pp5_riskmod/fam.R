# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Load functions for running models
source("../../_deps/run_model_functions.R")

# Load USC families
load("../usc_data/pp5/usc_families_mod.rData")

###############################################################################

# Split the families for parallelization
idx = split(1:length(usc_families_PP_mod), 
            ceiling(seq_along(1:length(usc_families_PP_mod)) / 
                      ceiling(length(usc_families_PP_mod)/20)))


# Run models for comparison on families
fam_output = lapply(idx[[s]], function(i){
  # Extract family
  fam_PP = usc_families_PP_mod[[i]]
  fam_BM_list = usc_families_BM_mod[[i]]
  
  # Run models
  out = try(run_models_5BC_rm(fam_PP, fam_BM_list), TRUE)
  return(out)
})


# Save results
save(fam_output, file = paste0("results/output/fam", s, ".rData"))
