library(PanelPRO)
library(BayesMendel)

# Load database that is back-compatible with the BayesMendel package
load("../../_deps/BackCompatibleDatabase.rData")

# Add PANC gene to the list of acceptable genes
NEW_GENE_TYPES = c(PanelPRO:::GENE_TYPES, "PANC")
NEW_GENE_TYPES = NEW_GENE_TYPES[!duplicated(NEW_GENE_TYPES)]
assignInNamespace("GENE_TYPES", NEW_GENE_TYPES, 
                  ns="PanelPRO", pos="package:PanelPRO")
NEW_ALL_GENE_VARIANT_TYPES = c(PanelPRO:::ALL_GENE_VARIANT_TYPES, 
                               PANC = "PANC")
NEW_ALL_GENE_VARIANT_TYPES = NEW_ALL_GENE_VARIANT_TYPES[!duplicated(NEW_ALL_GENE_VARIANT_TYPES)]
assignInNamespace("ALL_GENE_VARIANT_TYPES", NEW_ALL_GENE_VARIANT_TYPES, 
                  ns="PanelPRO", pos="package:PanelPRO")
NEW_DEFAULT_VARIANTS = c(PanelPRO:::DEFAULT_VARIANTS, 
                         PANC = "PANC")
NEW_DEFAULT_VARIANTS = NEW_DEFAULT_VARIANTS[!duplicated(NEW_DEFAULT_VARIANTS)]
assignInNamespace("DEFAULT_VARIANTS", NEW_DEFAULT_VARIANTS, 
                  ns="PanelPRO", pos="package:PanelPRO")


# Function to convert a family matrix from the PanelPro format to the 
# BayesMendel format
fam2BayesMendelFam = function(fam_PP) {
  
  fam_BM = fam_PP
  
  # Rename cancer affectedness columns
  names(fam_BM) = gsub("isAff", "Affected", names(fam_BM))
  
  # Rename cancer affected and age columns to match BayesMendel models
  names(fam_BM) = gsub("CBC", "BreastContralateral", names(fam_BM))
  names(fam_BM) = gsub("BC", "Breast", names(fam_BM))
  names(fam_BM) = gsub("COL", "Colon", names(fam_BM))
  names(fam_BM) = gsub("ENDO", "Endometrium", names(fam_BM))
  names(fam_BM) = gsub("MELA", "Skin", names(fam_BM))
  names(fam_BM) = gsub("OC", "Ovary", names(fam_BM))
  names(fam_BM) = gsub("PANC", "Pancreas", names(fam_BM))
  
  if (!("AgeBreastContralateral" %in% names(fam_BM))) {
    fam_BM$AgeBreastContralateral = 0
  }
  
  for (canc in c("Breast", "Ovary", "Colon", "Endometrium", "Pancreas", "Skin")) {
    if (paste0("Affected", canc) %in% names(fam_BM)) {
      fam_BM[fam_BM[paste0("Affected", canc)] == 0, paste0("Age", canc)] = 
        fam_BM[fam_BM[paste0("Affected", canc)] == 0, "CurAge"]
    }
  }
  
  # Rename the Sex column as Gender
  names(fam_BM) = gsub("Sex", "Gender", names(fam_BM))
  
  # Rename death columns (unused by BayesMendel models anyway)
  names(fam_BM) = gsub("isDead", "Death", names(fam_BM))
  names(fam_BM) = gsub("CurAge", "AgeDeath", names(fam_BM))
  
  # Set ethnic column
  fam_BM$ethnic = fam_PP$Ancestry
  
  # Need to set missing ages to 1 for MMRpro and melapro
  if ("AgeColon" %in% names(fam_BM)) {
    fam_BM$AgeColon[is.na(fam_BM$AgeColon)] = 1
    fam_BM$AgeEndometrium[is.na(fam_BM$AgeEndometrium)] = 1
  }
  if ("AgeSkin" %in% names(fam_BM)) {
    fam_BM$AgeSkin[is.na(fam_BM$AgeSkin)] = 1
  }
  
  return(fam_BM)
}


# Function to extract USC family structure from a single family
# Assumes that family has the following columns: 
# ID, MotherID, FatherID, Gender, isProband
# Returns a 4 x 2 matrix of 
# - the number of paternal aunts/uncles (patern.sibs)
# - the number of maternal aunts/uncles (matern.sibs)
# - the number of siblings of the proband (sibs)
# - the number children of the proband (ngchild)
# With respect to this language, the proband belongs to the "children" 
# generation of the family.
get_usc_fam_struct = function(usc_fam) {
  # Identify people who don't have parents
  isBothParentsMissing = (usc_fam$MotherID == 0) & (usc_fam$FatherID == 0)
  
  # Proband ID
  ID = usc_fam$ID[usc_fam$isProband==1]
  
  # IDs of mother and father
  mothID = usc_fam$MotherID[usc_fam$ID==ID]
  fathID = usc_fam$FatherID[usc_fam$ID==ID]
  
  # IDs of paternal grandparents
  mothID.patern = usc_fam$MotherID[usc_fam$ID==fathID]
  fathID.patern = usc_fam$FatherID[usc_fam$ID==fathID]
  
  # IDs of maternal grandparents
  mothID.matern = usc_fam$MotherID[usc_fam$ID==mothID]
  fathID.matern = usc_fam$FatherID[usc_fam$ID==mothID]
  
  # Number of paternal aunts and uncles
  nsibs.patern = c(sum((!(isBothParentsMissing) & usc_fam$MotherID == mothID.patern & 
                          usc_fam$FatherID == fathID.patern)[usc_fam$Sex==0], na.rm=TRUE), 
                   sum((!(isBothParentsMissing) & usc_fam$MotherID == mothID.patern & 
                          usc_fam$FatherID == fathID.patern)[usc_fam$Sex==1], na.rm=TRUE))
  # Number of maternal aunts and uncles
  nsibs.matern = c(sum((!(isBothParentsMissing) & usc_fam$MotherID == mothID.matern & 
                          usc_fam$FatherID == fathID.matern)[usc_fam$Sex==0], na.rm=TRUE), 
                   sum((!(isBothParentsMissing) & usc_fam$MotherID == mothID.matern & 
                          usc_fam$FatherID == fathID.matern)[usc_fam$Sex==1], na.rm=TRUE))
  # Number of paternal sisters and brothers
  nsibs = c(sum((!(isBothParentsMissing) & usc_fam$MotherID == mothID & 
                   usc_fam$FatherID == fathID)[usc_fam$Sex==0], na.rm=TRUE), 
            sum((!(isBothParentsMissing) & usc_fam$MotherID == mothID & 
                   usc_fam$FatherID == fathID)[usc_fam$Sex==1], na.rm=TRUE))
  # Number of grandchildren
  ngchild = c(sum((!(isBothParentsMissing) & usc_fam$MotherID == ID | 
                     usc_fam$FatherID == ID)[usc_fam$Sex==0], na.rm=TRUE), 
              sum((!(isBothParentsMissing) & usc_fam$MotherID == ID |
                     usc_fam$FatherID == ID)[usc_fam$Sex==1], na.rm=TRUE))
  
  # Drop the father, mother, and proband from counts, if necessary
  if (nsibs.patern[2] > 0) {
    nsibs.patern[2] = nsibs.patern[2]-1
  }
  if (nsibs.matern[1] > 0) {
    nsibs.matern[1] = nsibs.matern[1]-1
  }
  if (nsibs[usc_fam$Sex[usc_fam$ID==ID]+1] > 0) {
    nsibs[usc_fam$Sex[usc_fam$ID==ID]+1] = nsibs[usc_fam$Sex[usc_fam$ID==ID]+1]-1
  }
  
  # Return table of relative counts
  out = rbind(nsibs.patern, nsibs.matern, nsibs, ngchild)
  colnames(out) = c("female", "male")
  return(out)
}


# Extract posterior probabilities estimates as a named vector
# (ignores lower and upper CI)
extract_probs = function(res) {
  # Pull out posterior probabilities
  probs_df = res$posterior.prob[[1]]
  
  # Extract estimates only as a vector and name them using genotypes
  probs_vec = probs_df$estimate
  names(probs_vec) = probs_df$genes
  return(probs_vec)
}


# PanelPRO-7
my_PanelPRO7 = function(...) {
  PanelPRO::PanelPRO(genes=c("BRCA1", "BRCA2", "CDKN2A", 
                             "MLH1", "MSH2", "MSH6", "PANC"), 
                     cancers=c("Breast", "Ovarian", 
                               "Colorectal", "Endometrial", 
                               "Pancreas", "Melanoma"), 
                     ...)
}

# pancpro
my_pancpro = function(...) {
  PanelPRO::PanelPRO(genes="PANC", 
                     cancers="Pancreas", 
                     ...)
}

# melapro
my_melapro = function(...) {
  PanelPRO::PanelPRO(genes="CDKN2A", 
                     cancers="Melanoma", 
                     ...)
}



# Run PanelPRO-5BC and brcapro (using both PanelPRO and BayesMendel)
run_models_5BC = function(fam_PP, database = BackCompatibleDatabase) {
  
  # Convert PanelPRO-formatted family to BayesMendel format
  fam_BM = fam2BayesMendelFam(fam_PP)
  
  # Run PanelPRO-5BC
  probs_PP = extract_probs(PanelPRO::BRCAPRO5(fam_PP, database = database,
                                              max.mut = 2, parallel = FALSE, 
                                              allow.intervention = FALSE))
  # Run PanelPRO's brcapro
  probs_PP_brcapro = extract_probs(PanelPRO::BRCAPRO(fam_PP, database = database,
                                                     max.mut = 2, parallel = FALSE, 
                                                     allow.intervention = FALSE))
  
  # Run BayesMendel's brcapro
  dummy_comprisk = matrix(1, nrow=nrow(compriskSurv), ncol=ncol(compriskSurv))
  probs_BM_brcapro = as.numeric(BayesMendel::brcapro(fam_BM, 
                                                     fam_BM$ID[fam_BM$isProband==1], 
                                                     print = FALSE, 
                                                     params = brcaparams(
                                                       penetrance.net=penet.brca.net.PP,
                                                       comprisk=dummy_comprisk, 
                                                       CBCpenetrance.net=cbc.net.PP), 
                                                     imputeRelatives = FALSE)@posterior)
  names(probs_BM_brcapro) = names(probs_PP_brcapro)
  
  probs_PP_brcapro[setdiff(names(probs_PP), names(probs_PP_brcapro))] = NA
  probs_BM_brcapro[setdiff(names(probs_PP), names(probs_BM_brcapro))] = NA
  
  # Return posterior probabilities
  out = do.call(rbind, lapply(list(PanelPRO = probs_PP, 
                                   PP_brcapro = probs_PP_brcapro, 
                                   BM_brcapro = probs_BM_brcapro), 
                              function(x) x[match(names(probs_PP), names(x))]))
  
  # Rename genotypes
  colnames(out) = PanelPRO:::formatGeneNames(colnames(out), 
                                             format = "only_gene")
  return(out)
}


# Run PanelPRO-5BC and brcapro (using both PanelPRO and BayesMendel)
# Incorporate risk modifiers
run_models_5BC_rm = function(fam_PP, database = BackCompatibleDatabase) {
  
  # Convert PanelPRO-formatted family to BayesMendel format
  fam_BM = fam2BayesMendelFam(fam_PP)
  
  # Tumor biomarker testing results for brcapro
  brca.marker.testing = fam_PP[,c("ER", "CK14", "CK5.6", "PR", "HER2")]
  # Set 0 to 2 (negative test)
  brca.marker.testing[brca.marker.testing==0] = 2
  # Set NA to 0 (not tested)
  brca.marker.testing[is.na(brca.marker.testing)] = 0
  
  # Run PanelPRO-5BC
  probs_PP = extract_probs(PanelPRO::BRCAPRO5(fam_PP, database = database,
                                              max.mut = 2, parallel = FALSE, 
                                              allow.intervention = TRUE, 
                                              ignore.proband.germ = TRUE))
  # Run PanelPRO's brcapro
  probs_PP_brcapro = extract_probs(PanelPRO::BRCAPRO(fam_PP, database = database,
                                                     max.mut = 2, parallel = FALSE, 
                                                     allow.intervention = TRUE, 
                                                     ignore.proband.germ = TRUE))
  
  # Run BayesMendel's brcapro
  dummy_comprisk = matrix(1, nrow=nrow(compriskSurv), ncol=ncol(compriskSurv))
  probs_BM_brcapro = as.numeric(BayesMendel::brcapro(fam_BM, 
                                                     fam_BM$ID[fam_BM$isProband==1], 
                                                     marker.testing = brca.marker.testing,
                                                     print = FALSE, 
                                                     params = brcaparams(
                                                       penetrance.net=penet.brca.net.PP,
                                                       comprisk=dummy_comprisk, 
                                                       CBCpenetrance.net=cbc.net.PP), 
                                                     imputeRelatives=FALSE)@posterior)
  names(probs_BM_brcapro) = names(probs_PP_brcapro)
  
  probs_PP_brcapro[setdiff(names(probs_PP), names(probs_PP_brcapro))] = NA
  probs_BM_brcapro[setdiff(names(probs_PP), names(probs_BM_brcapro))] = NA
  
  # Return posterior probabilities
  out = do.call(rbind, lapply(list(PanelPRO = probs_PP, 
                                   PP_brcapro = probs_PP_brcapro, 
                                   BM_brcapro = probs_BM_brcapro), 
                              function(x) x[match(names(probs_PP), names(x))]))
  
  # Rename genotypes
  colnames(out) = PanelPRO:::formatGeneNames(colnames(out), 
                                             format = "only_gene")
  return(out)
}




# Run PanelPRO-7 and all sub-models (using both PanelPRO and BayesMendel)
run_models_7 = function(fam_PP, database = BackCompatibleDatabase) {
  
  # Convert PanelPRO-formatted family to BayesMendel format
  fam_BM = fam2BayesMendelFam(fam_PP)
  
  # Run PanelPRO-7
  probs_PP = extract_probs(my_PanelPRO7(fam_PP, database = database,
                                        max.mut = 2, parallel = FALSE, 
                                        allow.intervention = FALSE))
  # Run PanelPRO's brcapro
  probs_PP_brcapro = extract_probs(PanelPRO::BRCAPRO(fam_PP, database = database,
                                                     max.mut = 2, parallel = FALSE, 
                                                     allow.intervention = FALSE))
  
  # Run PanelPRO's mmrpro
  probs_PP_mmrpro = extract_probs(PanelPRO::MMRPRO(fam_PP, database = database,
                                                   max.mut = 2, parallel = FALSE, 
                                                   allow.intervention = FALSE))
  
  # Run PanelPRO's pancpro
  probs_PP_pancpro = extract_probs(my_pancpro(fam_PP, database = database,
                                              max.mut = 2, parallel = FALSE, 
                                              allow.intervention = FALSE))
  
  # Run PanelPRO's melapro
  probs_PP_melapro = extract_probs(my_melapro(fam_PP, database = database,
                                              max.mut = 2, parallel = FALSE, 
                                              allow.intervention = FALSE))
  
  
  # Run BayesMendel's brcapro
  dummy_comprisk = matrix(1, nrow=nrow(compriskSurv), ncol=ncol(compriskSurv))
  probs_BM_brcapro = as.numeric(BayesMendel::brcapro(fam_BM, 
                                                     fam_BM$ID[fam_BM$isProband==1], 
                                                     print = FALSE, 
                                                     params = brcaparams(
                                                       penetrance.net=penet.brca.net.PP,
                                                       comprisk=dummy_comprisk, 
                                                       CBCpenetrance.net=cbc.net.PP), 
                                                     imputeRelatives = FALSE)@posterior)
  names(probs_BM_brcapro) = names(probs_PP_brcapro)
  
  probs_PP_brcapro[setdiff(names(probs_PP), names(probs_PP_brcapro))] = NA
  probs_BM_brcapro[setdiff(names(probs_PP), names(probs_BM_brcapro))] = NA
  
  # Run BayesMendel's mmrpro
  probs_BM_mmrpro = as.numeric(BayesMendel::MMRpro(fam_BM, fam_BM$ID[fam_BM$isProband==1], 
                                                   print = FALSE, 
                                                   params = MMRparams(
                                                     penetrance.net=penet.mmr.net.PP),
                                                   imputeRelatives = FALSE)@posterior[1:2,1:2,1:2])[1:7]
  names(probs_BM_mmrpro) = names(probs_PP_mmrpro)[c(1,2,3,5,4,6,7)]
  
  probs_PP_mmrpro[setdiff(names(probs_PP), names(probs_PP_mmrpro))] = NA
  probs_BM_mmrpro[setdiff(names(probs_PP), names(probs_BM_mmrpro))] = NA
  
  # Run BayesMendel's pancpro
  probs_BM_pancpro = as.numeric(BayesMendel::pancpro(fam_BM, fam_BM$ID[fam_BM$isProband==1], 
                                                     params = pancparams(
                                                       penetrance.net=penet.panc.net.PP),
                                                     imputeRelatives = FALSE)@posterior[1:2])
  names(probs_BM_pancpro) = names(probs_PP_pancpro)
  
  probs_PP_pancpro[setdiff(names(probs_PP), names(probs_PP_pancpro))] = NA
  probs_BM_pancpro[setdiff(names(probs_PP), names(probs_BM_pancpro))] = NA
  
  # Run BayesMendel's melapro
  probs_BM_melapro = as.numeric(BayesMendel::melapro(fam_BM, fam_BM$ID[fam_BM$isProband==1], 
                                                     print = FALSE, 
                                                     params = melaparams(spm.lr.noncarrier=1, 
                                                                       spm.lr.carrier=1, 
                                                                       penetrance.net=penet.mela.hbi.net.PP), 
                                                     imputeRelatives=FALSE)@posterior[1:2])
  names(probs_BM_melapro) = names(probs_PP_melapro)
  
  probs_PP_melapro[setdiff(names(probs_PP), names(probs_PP_melapro))] = NA
  probs_BM_melapro[setdiff(names(probs_PP), names(probs_BM_melapro))] = NA
  
  # Return posterior probabilities
  out = do.call(rbind, lapply(list(PanelPRO = probs_PP, 
                                   PP_brcapro = probs_PP_brcapro, 
                                   PP_mmrpro = probs_PP_mmrpro, 
                                   PP_pancpro = probs_PP_pancpro, 
                                   PP_melapro = probs_PP_melapro, 
                                   BM_brcapro = probs_BM_brcapro, 
                                   BM_mmrpro = probs_BM_mmrpro, 
                                   BM_pancpro = probs_BM_pancpro, 
                                   BM_melapro = probs_BM_melapro), 
                              function(x) x[match(names(probs_PP), names(x))]))
  
  # Rename genotypes
  colnames(out) = PanelPRO:::formatGeneNames(colnames(out), 
                                             format = "only_gene")
  
  return(out)
}



# Run PanelPRO-11 and all sub-models (using both PanelPRO and BayesMendel)
run_models_11 = function(fam_PP, database = BackCompatibleDatabase) {
  
  # Convert PanelPRO-formatted family to BayesMendel format
  fam_BM = fam2BayesMendelFam(fam_PP)
  
  # Run PanelPRO-11
  probs_PP = extract_probs(PanelPRO::PanelPRO11(fam_PP, database = database,
                                                max.mut = 2, parallel = FALSE, 
                                                allow.intervention = FALSE))
  # Run PanelPRO's brcapro
  probs_PP_brcapro = extract_probs(PanelPRO::BRCAPRO(fam_PP, database = database,
                                                     max.mut = 2, parallel = FALSE, 
                                                     allow.intervention = FALSE))
  
  # Run PanelPRO's mmrpro
  probs_PP_mmrpro = extract_probs(PanelPRO::MMRPRO(fam_PP, database = database,
                                                   max.mut = 2, parallel = FALSE, 
                                                   allow.intervention = FALSE))
  
  # Run PanelPRO's melapro
  probs_PP_melapro = extract_probs(my_melapro(fam_PP, database = database,
                                              max.mut = 2, parallel = FALSE, 
                                              allow.intervention = FALSE))
  
  
  # Run BayesMendel's brcapro
  dummy_comprisk = matrix(1, nrow=nrow(compriskSurv), ncol=ncol(compriskSurv))
  probs_BM_brcapro = as.numeric(BayesMendel::brcapro(fam_BM, 
                                                     fam_BM$ID[fam_BM$isProband==1], 
                                                     print = FALSE, 
                                                     params = brcaparams(
                                                       penetrance.net=penet.brca.net.PP,
                                                       comprisk=dummy_comprisk, 
                                                       CBCpenetrance.net=cbc.net.PP), 
                                                     imputeRelatives = FALSE)@posterior)
  names(probs_BM_brcapro) = names(probs_PP_brcapro)
  
  probs_PP_brcapro[setdiff(names(probs_PP), names(probs_PP_brcapro))] = NA
  probs_BM_brcapro[setdiff(names(probs_PP), names(probs_BM_brcapro))] = NA
  
  # Run BayesMendel's mmrpro
  probs_BM_mmrpro = as.numeric(BayesMendel::MMRpro(fam_BM, fam_BM$ID[fam_BM$isProband==1], 
                                                   print = FALSE, 
                                                   params = MMRparams(
                                                     penetrance.net=penet.mmr.net.PP),
                                                   imputeRelatives = FALSE)@posterior[1:2,1:2,1:2])[1:7]
  names(probs_BM_mmrpro) = names(probs_PP_mmrpro)[c(1,2,3,5,4,6,7)]
  
  probs_PP_mmrpro[setdiff(names(probs_PP), names(probs_PP_mmrpro))] = NA
  probs_BM_mmrpro[setdiff(names(probs_PP), names(probs_BM_mmrpro))] = NA
  
  # Run BayesMendel's melapro
  probs_BM_melapro = as.numeric(BayesMendel::melapro(fam_BM, fam_BM$ID[fam_BM$isProband==1], 
                                                     print = FALSE, 
                                                     params = melaparams(spm.lr.noncarrier=1, 
                                                                       spm.lr.carrier=1, 
                                                                       penetrance.net=penet.mela.hbi.net.PP), 
                                                     imputeRelatives = FALSE)@posterior[1:2])
  names(probs_BM_melapro) = names(probs_PP_melapro)
  
  probs_PP_melapro[setdiff(names(probs_PP), names(probs_PP_melapro))] = NA
  probs_BM_melapro[setdiff(names(probs_PP), names(probs_BM_melapro))] = NA
  
  # Return posterior probabilities
  out = do.call(rbind, lapply(list(PanelPRO = probs_PP, 
                                   PP_brcapro = probs_PP_brcapro, 
                                   PP_mmrpro = probs_PP_mmrpro, 
                                   PP_melapro = probs_PP_melapro, 
                                   BM_brcapro = probs_BM_brcapro, 
                                   BM_mmrpro = probs_BM_mmrpro, 
                                   BM_melapro = probs_BM_melapro), 
                              function(x) x[match(names(probs_PP), names(x))]))
  
  # Rename genotypes
  colnames(out) = PanelPRO:::formatGeneNames(colnames(out), 
                                             format = "only_gene")
  
  return(out)
}


# Run PanelPRO-11 and all sub-models (using both PanelPRO and BayesMendel)
# Incorporate risk modifiers
run_models_11_rm = function(fam_PP, database = BackCompatibleDatabase) {
  
  # Convert PanelPRO-formatted family to BayesMendel format
  fam_BM = fam2BayesMendelFam(fam_PP)
  
  # Tumor biomarker testing results for brcapro
  brca.marker.testing = fam_PP[,c("ER", "CK14", "CK5.6", "PR", "HER2")]
  # Set 0 to 2 (negative test)
  brca.marker.testing[brca.marker.testing==0] = 2
  # Set NA to 0 (not tested)
  brca.marker.testing[is.na(brca.marker.testing)] = 0
  
  # Tumor biomarker testing results for MMRpro
  mmr.marker.testing = data.frame(MSI=fam_PP$MSI, location=0)
  # Set 0 to 2 (negative test)
  mmr.marker.testing$MSI[mmr.marker.testing$MSI==0] = 2
  # Set NA to 0 (not tested)
  mmr.marker.testing$MSI[is.na(mmr.marker.testing$MSI)] = 0
  
  # Run PanelPRO-11
  probs_PP = extract_probs(PanelPRO::PanelPRO11(fam_PP, database = database,
                                                max.mut = 2, parallel = FALSE, 
                                                allow.intervention = TRUE, 
                                                ignore.proband.germ = TRUE))
  # Run PanelPRO's brcapro
  probs_PP_brcapro = extract_probs(PanelPRO::BRCAPRO(fam_PP, database = database,
                                                     max.mut = 2, parallel = FALSE, 
                                                     allow.intervention=TRUE, 
                                                     ignore.proband.germ = TRUE))
  
  # Run PanelPRO's mmrpro
  probs_PP_mmrpro = extract_probs(PanelPRO::MMRPRO(fam_PP, database = database,
                                                   max.mut = 2, parallel = FALSE, 
                                                   allow.intervention = TRUE, 
                                                   ignore.proband.germ = TRUE))
  
  # Run PanelPRO's melapro
  probs_PP_melapro = extract_probs(my_melapro(fam_PP, database = database,
                                              max.mut = 2, parallel = FALSE, 
                                              allow.intervention = TRUE, 
                                              ignore.proband.germ = TRUE))
  
  
  # Run BayesMendel's brcapro
  dummy_comprisk = matrix(1, nrow=nrow(compriskSurv), ncol=ncol(compriskSurv))
  probs_BM_brcapro = as.numeric(BayesMendel::brcapro(fam_BM, 
                                                     fam_BM$ID[fam_BM$isProband==1], 
                                                     marker.testing = brca.marker.testing,
                                                     print = FALSE, 
                                                     params = brcaparams(
                                                       penetrance.net=penet.brca.net.PP,
                                                       comprisk=dummy_comprisk, 
                                                       CBCpenetrance.net=cbc.net.PP), 
                                                     imputeRelatives = FALSE)@posterior)
  names(probs_BM_brcapro) = names(probs_PP_brcapro)
  
  probs_PP_brcapro[setdiff(names(probs_PP), names(probs_PP_brcapro))] = NA
  probs_BM_brcapro[setdiff(names(probs_PP), names(probs_BM_brcapro))] = NA
  
  # Run BayesMendel's mmrpro
  probs_BM_mmrpro = as.numeric(BayesMendel::MMRpro(fam_BM, fam_BM$ID[fam_BM$isProband==1], 
                                                   marker.testing = mmr.marker.testing,
                                                   print = FALSE, 
                                                   params = MMRparams(
                                                     penetrance.net=penet.mmr.net.PP),
                                                   imputeRelatives = FALSE)@posterior[1:2,1:2,1:2])[1:7]
  names(probs_BM_mmrpro) = names(probs_PP_mmrpro)[c(1,2,3,5,4,6,7)]
  
  probs_PP_mmrpro[setdiff(names(probs_PP), names(probs_PP_mmrpro))] = NA
  probs_BM_mmrpro[setdiff(names(probs_PP), names(probs_BM_mmrpro))] = NA
  
  # Run BayesMendel's melapro
  probs_BM_melapro = as.numeric(BayesMendel::melapro(fam_BM, fam_BM$ID[fam_BM$isProband==1], 
                                                     print = FALSE, 
                                                     params = melaparams(spm.lr.noncarrier=1, 
                                                                       spm.lr.carrier=1, 
                                                                       penetrance.net=penet.mela.hbi.net.PP), 
                                                     imputeRelatives = FALSE)@posterior[1:2])
  names(probs_BM_melapro) = names(probs_PP_melapro)
  
  probs_PP_melapro[setdiff(names(probs_PP), names(probs_PP_melapro))] = NA
  probs_BM_melapro[setdiff(names(probs_PP), names(probs_BM_melapro))] = NA
  
  # Return posterior probabilities
  out = do.call(rbind, lapply(list(PanelPRO = probs_PP, 
                                   PP_brcapro = probs_PP_brcapro, 
                                   PP_mmrpro = probs_PP_mmrpro, 
                                   PP_melapro = probs_PP_melapro, 
                                   BM_brcapro = probs_BM_brcapro, 
                                   BM_mmrpro = probs_BM_mmrpro, 
                                   BM_melapro = probs_BM_melapro), 
                              function(x) x[match(names(probs_PP), names(x))]))
  
  # Rename genotypes
  colnames(out) = PanelPRO:::formatGeneNames(colnames(out), 
                                             format = "only_gene")
  
  return(out)
}
