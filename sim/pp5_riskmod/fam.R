# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Load functions for running models
source("../../_deps/run_model_functions.R")

# Load simulation functions
files.sources = paste0("../../simulate_families/", 
                       list.files("../../simulate_families/"))
sapply(files.sources, source)

# Load USC families
load("../../usc/usc_data/pp5/usc_families.rData")

###############################################################################

## Cancers and genes to use
# Short cancer names (including CBC)
cancers = c("BC", "OC", "CBC")
# Look up long cancer names (don't include CBC)
cancers_long = PanelPRO:::CANCER_NAME_MAP$long[sapply(cancers[-length(cancers)], function(x){
  which(x==PanelPRO:::CANCER_NAME_MAP$short)
})]
# Genes
genes = c("ATM", "BRCA1", "BRCA2", "CHEK2", "PALB2")


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
                         ppd=BackCompatibleDatabase)
dummy.db$Contralateral = BackCompatibleDatabase$Contralateral

# Run `checkFam` on the dummy family
dummy.fam.checked = checkFam(dummy.fam, dummy.db)$ped_list[[1]]

# Cancer penetrance densities and survivals
CP = calcCancerPenetrance(dummy.fam.checked, dummy.db, 
                          max_mut=2, net=TRUE, consider.modification=FALSE)

# Extract allele frequencies from database
alleleFreq = BackCompatibleDatabase$AlleleFrequency[,"nonAJ"][genes]

###############################################################################

# 1000 family simulations
nsim = 1000

# Run models for comparison on simulated families
fam_output = lapply(1:nsim, function(i){
  # Simulate family
  
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
  fam = try(sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                       alleleFreq, CP, genes, cancers_long, 
                       includeGeno=TRUE, 
                       BiomarkerTesting=PanelPRODatabase$BiomarkerTesting, 
                       includeBiomarkers=TRUE), TRUE)
  if (class(fam)=="try-error") {
    fam = try(sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                         alleleFreq, CP, genes, cancers_long, 
                         includeGeno=TRUE, 
                         BiomarkerTesting=PanelPRODatabase$BiomarkerTesting, 
                         includeBiomarkers=TRUE))
  }
  
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
  # Put together list with modifiers
  fam_BM_list = list(fam = fam_BM,
                     brcapro = list(germline.testing = NULL,
                                    marker.testing = brca.marker.testing))
  
  # Run models
  out = try(run_models_5BC_rm(fam_PP, fam_BM_list), TRUE)
  return(list(fam = fam, probs = out))
})


# Save results
save(fam_output, file = paste0("results/output/fam", s, ".rData"))
