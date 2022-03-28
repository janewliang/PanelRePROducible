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
# 40 different loci, each with MAF 0.1 (so prob of carrying = 0.19)
# "latent score" is the arithmetic sum, so it ranges from 0 to 40

# Expected proportion of each score
# score_dist = sapply(0:40, function(k) { choose(40, k)*(0.19^k)*(0.81^(40-k)) })
# barplot(score_dist, names.arg = 0:40)
# sum(score_dist[1:8]) # median occurs right between scores 7 and 8

# Additional risk factor (multiply genotype-specific cancer penetrance by): 
# 0 (0.0002) 0.6
# 1 (0.0020) 0.6 
# 2 (0.0094) 0.6
# 3 (0.0279) 0.65 
# 4 (0.0604) 0.7
# 5 (0.1021) 0.75 
# 6 (0.1397) 0.85
# 7 (0.1592) 0.95 
# 8 (0.1540) 1.05
# 9 (0.1284) 1.15
# 10 (0.0934) 1.25
# 11 (0.0597) 1.3
# 12 (0.0339) 1.35
# 13+ (0.02968973) 1.4 
# Divide all factors by 1.004991 so that weighted mean close to 1
# (i.e. total number of cancer cases should stay about the same)
latent = data.frame(af = 0.1, 
                    fact = c(rep(0.6, 3), 0.65, 0.7, 0.75, 0.85, 0.95, 
                             1.05, 1.15, 1.25, 1.3, 1.35, rep(1.4, 28)) / 1.004991)

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
