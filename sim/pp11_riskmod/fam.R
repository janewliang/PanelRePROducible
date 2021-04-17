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
  fam = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                      BackCompatibleDatabase, genes, cancers, 
                      includeGeno=TRUE, includeBiomarkers=TRUE)
  
  # PanelPRO family
  fam_PP = fam
  # Remove mutations
  fam_PP[,genes] = NULL
  
  # BayesMendel family
  fam_BM = fam2BayesMendelFam(fam)
  # brcapro tumor marker testing
  brca.marker.testing = fam[c("ER", "CK14", "CK5.6", "PR", "HER2")]
  # Set 0 to 2 (negative test)
  brca.marker.testing[brca.marker.testing==0] = 2
  # Set NA to 0 (not tested)
  brca.marker.testing[is.na(brca.marker.testing)] = 0
  if (all(rowSums(brca.marker.testing)==0)) {
    brca.marker.testing = NULL
  }
  # mmrpro tumor marker testing
  mmr.marker.testing = data.frame(MSI=fam$MSI, location=0)
  # Set 0 to 2 (negative test)
  mmr.marker.testing[mmr.marker.testing==0] = 2
  # Set NA to 0 (not tested)
  mmr.marker.testing[is.na(mmr.marker.testing)] = 0
  if (all(rowSums(mmr.marker.testing)==0)) {
    mmr.marker.testing = NULL
  }
  # Put together list with modifiers
  fam_BM_list = list(fam = fam_BM,
                     brcapro=list(germline.testing=NULL, 
                                  marker.testing=brca.marker.testing), 
                     mmrpro=list(germline.testing=NULL, 
                                 marker.testing=mmr.marker.testing), 
                     melapro=list(germline.testing=NULL))
  
  # Run models
  out = run_models_11_rm(fam_PP, fam_BM_list)
  return(list(fam = fam, probs = out))
})


# Save results
save(fam_output, file = paste0("results/output/fam", s, ".rData"))
