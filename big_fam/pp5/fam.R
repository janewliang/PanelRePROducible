# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Load functions for running models
source("../../_deps/run_model_functions.R")

# Load simulation functions
files.sources = paste0("../../simulate_families/", 
                       list.files("../../simulate_families/"))
sapply(files.sources, source)

###############################################################################

## Cancers and genes to use
# Cancers
cancers = c("Breast", "Ovarian")
# Genes
genes = c("ATM", "BRCA1", "BRCA2", "CHEK2", "PALB2")

###############################################################################

# 1000 family simulations
nsim = 1000

# Run models for comparison on simulated families
fam_output = lapply(1:nsim, function(i){
  
  # Set random seed
  set.seed(s*nsim+i)
  
  # Assign number of males and females in each branch of the family to 4.
  # Paternal aunt, paternal uncles
  nSibsPatern = c(4, 4)
  # Maternal aunts, maternal uncles
  nSibsMatern = c(4, 4)
  # Sisters and brothers
  nSibs = c(4, 4)
  # We make the assumption that the number of sons and daughters for the 
  # proband and all siblings, is the same. Nieces and nephews of the proband 
  # are not sampled separately.
  nGrandchild = c(4, 4)
  
  
  # Simulate family
  fam_PP = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                         BackCompatibleDatabase, genes, cancers, 
                         includeGeno=TRUE)
  fam_BM = fam2BayesMendelFam(fam_PP)
  
  # Run models
  out = run_models_5BC(fam_PP, fam_BM)
  return(list(fam = fam_PP, probs = out))
})


# Save results
save(fam_output, file = paste0("results/output/fam", s, ".rData"))
