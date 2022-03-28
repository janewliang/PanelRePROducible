# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Load functions for running models
source("../../_deps/run_model_functions.R")

# Load simulation functions
files.sources = paste0("../../simulate_families_latent/", 
                       list.files("../../simulate_families_latent/"))
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

# Set up data frame for latent risk factor
# 20 different loci, each with MAF 0.1 (so prob of carrying = 0.19)
# "latent score" is the arithmetic mean, so it ranges from 0 to 20
# Additional risk factor (multiply genotype-specific cancer penetrance by): 
# 0 (0.01478088) 0.8 * 1.010572
# 1 (0.06934241) 0.85 * 1.010572
# 2 (0.1545223) 0.9 * 1.010572
# 3 (0.2174758) 0.95 * 1.010572
# 4 (0.2168046) 1 * 1.010572
# 5 (0.1627373) 1.05 * 1.010572
# 6 (0.09543235) 1.1 * 1.010572
# 7 (0.04477073) 1.15 * 1.010572
# 8+ (0.02413363) 1.2 * 1.010572
# weighted mean is 1 (i.e. total number of cancer cases should stay about the same)
latent = data.frame(af = 0.1, 
                    fact = c(0.80, 0.85, 0.90, 0.95, 1.00, 
                             1.05, 1.10, 1.15, rep(1.2, 13)) * 1.010572)

###############################################################################

# 1000 family simulations
nsim = 1000

# Run models for comparison on simulated families
fam_output = lapply(1:nsim, function(i){
  
  # Set random seed
  set.seed(s*nsim+i)
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
                         includeGeno = TRUE, ageMin = 1, latent = latent)
  
  # Run models
  out = run_models_11(fam_PP)
  return(list(fam = fam_PP, probs = out))
})

# Save results
save(fam_output, file = paste0("results/output/fam", s, ".rData"))
