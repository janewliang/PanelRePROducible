# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Load functions for running models
source("../../_deps/run_model_functions.R")

# Load simulation functions
files.sources = paste0("../../simulate_families/", 
                       list.files("../../simulate_families/"))
sapply(files.sources, source)

# Load USC families
load("../../usc/usc_data/pp7/usc_families.rData")

###############################################################################

# Create a new database that is similar to BackCompatibleDatabase, but does
# not take into account penetrance information that is not contained in the
# BayesMendel models (e.g. pancreatic cancer penetrance for BRCA1 carriers). 
# This database will be used to simulate the families. 

# Duplicate database
BackCompatibleDatabase_sub = BackCompatibleDatabase

# Replace penetrances with SEER unless they belong in the BayesMendel models
for (cancer in dimnames(BackCompatibleDatabase_sub$Penetrance)$Cancer) {
  for (sex in c("Female", "Male")) {
    if (cancer == "Breast" || cancer == "Ovarian") {
      keep_genes = c("BRCA1", "BRCA2")
    } else if (cancer == "Colorectal" || cancer == "Endometrial") {
      keep_genes = c("MLH1", "MSH2", "MSH6")
    } else if (cancer == "Pancreas") {
      keep_genes = c("PANC")
    } else if (cancer == "Melanoma") {
      keep_genes = c("CDKN2A")
    } else {
      keep_genes = c()
    }
    is_overwrite = !(PanelPRO:::formatGeneNames(
      dimnames(BackCompatibleDatabase_sub$Penetrance)$Gene, 
      format = "only_gene") %in% keep_genes)
    if (any(is_overwrite == TRUE)) {
      BackCompatibleDatabase_sub$Penetrance[cancer,is_overwrite,"All_Races",sex,,"Net"] =
        matrix(BackCompatibleDatabase_sub$Penetrance[cancer,"SEER","All_Races",sex,,"Net"],
               nrow=sum(is_overwrite), 
               ncol=dim(BackCompatibleDatabase_sub$Penetrance)[5], 
               byrow=TRUE)
    }
  }
}

###############################################################################

## Cancers and genes to use
# Cancers
cancers = c("Breast", "Colorectal", "Endometrial", 
            "Melanoma", "Ovarian", "Pancreas")
# Genes
genes = c("BRCA1", "BRCA2", "CDKN2A", "MLH1", "MSH2", "MSH6", "PANC")

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
                         BackCompatibleDatabase_sub, genes, cancers, 
                         includeGeno=TRUE)
  
  # Run models
  out = run_models_7(fam_PP)
  return(list(fam = fam_PP, probs = out))
})


# Save results
save(fam_output, file = paste0("results/output/fam", s, ".rData"))
