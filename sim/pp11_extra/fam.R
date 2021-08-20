# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Load functions for running models
source("../../_deps/run_model_functions.R")

# Load simulation functions
files.sources = paste0("../../simulate_families/", 
                       list.files("../../simulate_families/"))
sapply(files.sources, source)

# Load USC families
load("../../usc/usc_data/pp11/usc_families.rData")

###############################################################################

## Cancers and genes to use
# Cancers
cancers = c("Brain", "Breast", "Colorectal", "Endometrial", 
            "Gastric", "Kidney", "Melanoma", "Ovarian", 
            "Pancreas", "Prostate", "Small Intestine")
# Genes
genes = c("ATM", "BRCA1", "BRCA2", "CDKN2A", "CHEK2", "EPCAM", 
          "MLH1", "MSH2", "MSH6", "PALB2", "PMS2")

###############################################################################

# 1000 family simulations
nsim = 1000

# Run models for comparison on simulated families
fam_output = lapply(1:nsim, function(i){
  
  # Set random seed
  set.seed((1000+s)*nsim+i)
  # Sample from USC families
  rn = sample(1:length(usc_families_PP), 1)
  
  # Get USC family structure
  famStruct = get_usc_fam_struct(usc_families_PP[[rn]])
  
  # Randomly sample number of males and females in each branch of the family.
  # Paternal aunt, paternal uncles
  nSibsPatern = famStruct["nsibs.patern",]
  # Maternal aunts, maternal uncles
  nSibsMatern = famStruct["nsibs.matern",]
  # Sisters and brothers
  nSibs = famStruct["nsibs",]
  # We make the assumption that the number of sons and daughters for the 
  # proband and all siblings, is the same. Nieces and nephews of the proband 
  # are not sampled separately.
  nGrandchild = famStruct["ngchild",]
  
  
  # Simulate family
  fam_PP = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                         BackCompatibleDatabase, genes, cancers, 
                         includeGeno=TRUE)
  
  # Run models
  out = run_models_11(fam_PP)
  return(list(fam = fam_PP, probs = out))
})


# Save results
save(fam_output, file = paste0("results/output/fam", s, ".rData"))
