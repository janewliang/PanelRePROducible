# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Load functions for running models
source("../../_deps/run_model_functions.R")

# Load simulation functions
files.sources = paste0("../../_deps/family_simulation_functions/", 
                       list.files("../../_deps/family_simulation_functions/"))
sapply(files.sources, source)

# Load HCP families
load("../../hcp/hcp_data/pp7/hcp_families.rData")

###############################################################################

# Create a new database that is similar to BackCompatibleDatabase, but does
# not take into account penetrance information that is not contained in the
# BayesMendel models (e.g. pancreatic cancer penetrance for BRCA1 carriers)
# This database should be used to simulate families

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
    is_overwite = !(dimnames(BackCompatibleDatabase_sub$Penetrance)$Gene %in% keep_genes)
    if (any(is_overwite == TRUE)) {
      BackCompatibleDatabase_sub$Penetrance[cancer,is_overwite,"All_Races",sex,,"Net"] =
        matrix(BackCompatibleDatabase_sub$Penetrance[cancer,"SEER","All_Races",sex,,"Net"],
               nrow=sum(is_overwite), 
               ncol=dim(BackCompatibleDatabase_sub$Penetrance)[5], 
               byrow=TRUE)
    }
  }
}

###############################################################################

## Cancers and genes to use
# Short cancer names (including CBC)
cancers = c("BC", "COL", "ENDO", "MELA", "OC", "PANC", "CBC")
# Look up long cancer names (don't include CBC)
cancers_long = PanelPRO:::CANCER_NAME_MAP$long[sapply(cancers[-length(cancers)], function(x){
  which(x==PanelPRO:::CANCER_NAME_MAP$short)
})]
# Genes
genes = c("BRCA1", "BRCA2", "CDKN2A", "MLH1", "MSH2", "MSH6", "PANC")


## Create a dummy family to generate penetrance densities and survivals
# Empty data frame of cancer affectation statuses
dummy.cancers = setNames(as.data.frame(matrix(0, nrow=2, ncol=length(cancers))), 
                         paste0("isAff", cancers))
# Empty data frame of cancer ages
dummy.ages = setNames(as.data.frame(matrix(NA, nrow=2, ncol=length(cancers))), 
                      paste0("Age", cancers))
# Dummy family
dummy.fam = data.frame(ID=c(1,2), 
                       MotherID=c(0,1), 
                       FatherID=c(0,0), 
                       Sex=c(0,1), 
                       isProband=c(0,1), 
                       Twins=c(0,0), 
                       Ancestry=rep("nonAJ", 2), 
                       CurAge=c(60,30), 
                       isDead=c(0,0), 
                       race=rep("All_Races", 2), 
                       dummy.cancers, 
                       dummy.ages)
# Assign no risk modifiers
dummy.fam$interAge = dummy.fam$riskmod = list(character(0))


## Get cancer penetrance densities and survivals from the dummy family
# Build a dummy database
dummy.db = buildDatabase(genes=genes, 
                         cancers=cancers_long, 
                         ppd=BackCompatibleDatabase_sub)
dummy.db$Contralateral = BackCompatibleDatabase_sub$Contralateral

# Run `checkFam` on the dummy family
dummy.fam.checked = checkFam(dummy.fam, dummy.db)$ped_list[[1]]

# Cancer penetrance densities and survivals
CP = calcCancerPenetrance(dummy.fam.checked, dummy.db, 
                          max_mut=2, net=TRUE, consider.modification=FALSE)

# Extract allele frequencies from database
alleleFreq = BackCompatibleDatabase_sub$AlleleFrequency[,"nonAJ"][genes]

###############################################################################

# 1000 family simulations
nsim = 1000

# Run models for comparison on simulated families
fam_output = lapply(1:nsim, function(i){
  # Simulate family
  
  # Set random seed
  set.seed(s*nsim+i)
  # Sample from HCP families
  rn = sample(1:length(hcp_families_PP), 1)
  
  # Get HCP family structure
  famStruct = get_hcp_fam_struct(hcp_families_PP[[rn]])
  
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
  fam_PP = try(sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                          alleleFreq, CP, genes, cancers_long, 
                          includeGeno=TRUE), TRUE)
  if (class(fam_PP)=="try-error") {
    fam_PP = try(sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                            alleleFreq, CP, genes, cancers_long, 
                            includeGeno=TRUE))
  }
  fam_BM = fam2BayesMendelFam(fam_PP)
  
  # Run models
  out = try(run_models_7(fam_PP, fam_BM), TRUE)
  return(list(fam = fam_PP, probs = out))
})


# Save results
save(fam_output, file = paste0("results/output/fam", s, ".rData"))
