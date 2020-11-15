# Load HCP data
load("../hcp_data.RData")

# 11 mutations to consider
mutations = c("ATM", "BRCA1", "BRCA2", "CDKN2A", "CHEK2", "EPCAM",
              "MLH1", "MSH2", "MSH6", "PALB2", "PMS2")

# Drop probands who are carriers of other mutations, but not one of the 11
isCarrOtherMut = apply((hcp[,paste0("res", 1:5)] == 1) & 
                         !(apply(hcp[,paste0("gene", 1:5)], 1:2, function(x){
                           x %in% c(mutations, NA)})), 
                       1, any) &
  !(apply((hcp[,paste0("res", 1:5)] == 1) & 
            apply(hcp[,paste0("gene", 1:5)], 1:2, function(x){
              x %in% mutations}), 
          1, any))
isCarrOtherMut[is.na(isCarrOtherMut)==TRUE] = FALSE
hcp = hcp[!(hcp$Proband==1 & isCarrOtherMut==TRUE),]

# Carriers = at least one pathogenic mutation
# Noncarriers = negative for all mutations
# Drop probands with VUSs and no pathogenic genes
isVUS = 
  # Individual has VUSs
  apply(hcp[,paste0("res", 1:5)], 1, function(x){ 0 %in% x}) & 
  # Individual has no pathogenic genes
  !apply((hcp[,paste0("res", 1:5)] == 1) & # Pathogenic
          apply(hcp[,paste0("gene", 1:5)], 1:2, function(x){ # One of the genes
            x %in% mutations}), 
        1, any)
hcp = hcp[!(hcp$Proband==1 & isVUS==TRUE),]

# Create data frame of mutations
hcp_muts = data.frame(sapply(mutations, function(mut){
  (hcp$gene1==mut & hcp$res1==1) | (hcp$gene2==mut & hcp$res2==1) | 
    (hcp$gene3==mut & hcp$res3==1) | (hcp$gene4==mut & hcp$res4==1) | 
    (hcp$gene5==mut & hcp$res5==1)
}))
hcp_muts[is.na(hcp_muts)] = 0

# 11 cancers to consider
# Create data frame of cancers
hcp_cancers = hcp[,c("AffectedBrain", "AgeBrain", 
                     "AffectedBreast", "AgeBreast", 
                     "AffectedColon", "AgeColon", 
                     "AffectedEndometrium", "AgeEndometrium", 
                     "AffectedGastric", "AgeGastric", 
                     "AffectedKidney", "AgeKidney", 
                     "AffectedSkin", "AgeSkin", 
                     "AffectedOvary", "AgeOvary", 
                     "AffectedPancreas", "AgePancreas", 
                     "AffectedProstate", "AgeProstate", 
                     "AffectedSmallIntestine", "AgeSmallIntestine")]


# Markers to consider 
markers = c("ER", "CK14", "CK5.6", "PR", "HER2", "MSI")
# Tumor markers for BayesMendel 
hcp_markers_BM = hcp[markers]
# Tumor markers for PanelPRO 
hcp_markers_PP = hcp_markers_BM
# Set 0 to NA (not tested)
hcp_markers_PP[hcp_markers_PP==0] = NA
# Set 2 to 0 (negative test)
hcp_markers_PP[hcp_markers_PP==2] = 0


# Build data frame of all families
hcp_df = data.frame(FamID=hcp$FamID, 
                    ID = hcp$ID, 
                    MotherID = hcp$MotherID, 
                    FatherID = hcp$FatherID, 
                    Gender = hcp$Gender, 
                    isProband = hcp$Proband, 
                    Twins = 0, 
                    ethnic = "nonAJ",
                    AgeDeath = as.numeric(hcp$Age.of.entry), 
                    hcp_cancers, 
                    Death = as.numeric(hcp$Deceased.Status), 
                    hcp_muts, stringsAsFactors=FALSE)
# Identify individuals with missing values
drop_idx = rowSums(is.na(hcp_df)) > 0
# Date for germline testing
germline.date = as.Date(hcp$GenTestTable.Result.Date, format="%Y-%m-%d")

# Set maximum age to PanelPRO's MAXAGE (94)
hcp_df$AgeDeath[hcp_df$AgeDeath > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeBrain[hcp_df$AgeBrain > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeBreast[hcp_df$AgeBreast > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeColon[hcp_df$AgeColon > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeEndometrium[hcp_df$AgeEndometrium > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeGastric[hcp_df$AgeGastric > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeKidney[hcp_df$AgeKidney > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeSkin[hcp_df$AgeSkin > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeOvary[hcp_df$AgeOvary > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgePancreas[hcp_df$AgePancreas > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeProstate[hcp_df$AgeProstate > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeSmallIntestine[hcp_df$AgeSmallIntestine > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
# Hack to substitute affected age=1 with 2
# BM models use 1 to code for unknown affected age
hcp_df$AgeBreast[hcp_df$AgeBreast==1 & hcp_df$AffectedBreast==1] = 2
hcp_df$AgeOvary[hcp_df$AgeOvary==1 & hcp_df$AffectedOvary==1] = 2
hcp_df$AgeColon[hcp_df$AgeColon==1 & hcp_df$AffectedColon==1] = 2
hcp_df$AgeEndometrium[hcp_df$AgeEndometrium==1 & hcp_df$AffectedEndometrium==1] = 2
hcp_df$AgeSkin[hcp_df$AgeSkin==1 & hcp_df$AffectedSkin==1] = 2


# BayesMendel without modifiers
hcp_BM_df = data.frame(hcp_df, AgeBreastContralateral=0, 
                       germline.date=germline.date)
# BayesMendel with modifiers
hcp_BM_mod_df = data.frame(hcp_BM_df, hcp_markers_BM)

# PanelPRO
hcp_PP_df = data.frame(hcp_df, race="All_Races", 
                       germline.date=germline.date, 
                       stringsAsFactors=FALSE)
# Re-name columns
names(hcp_PP_df) = gsub("Gender", "Sex", names(hcp_PP_df))
names(hcp_PP_df) = gsub("AgeDeath", "CurAge", names(hcp_PP_df))
names(hcp_PP_df) = gsub("Death", "isDead", names(hcp_PP_df))
names(hcp_PP_df) = gsub("Affected", "isAff", names(hcp_PP_df))
names(hcp_PP_df) = gsub("ethnic", "Ancestry", names(hcp_PP_df))
# Use cancer codes
names(hcp_PP_df) = gsub("Brain", "BRA", names(hcp_PP_df))
names(hcp_PP_df) = gsub("Breast", "BC", names(hcp_PP_df))
names(hcp_PP_df) = gsub("Colon", "COL", names(hcp_PP_df))
names(hcp_PP_df) = gsub("Endometrium", "ENDO", names(hcp_PP_df))
names(hcp_PP_df) = gsub("Gastric", "GAS", names(hcp_PP_df))
names(hcp_PP_df) = gsub("Kidney", "KID", names(hcp_PP_df))
names(hcp_PP_df) = gsub("Skin", "MELA", names(hcp_PP_df))
names(hcp_PP_df) = gsub("Ovary", "OC", names(hcp_PP_df))
names(hcp_PP_df) = gsub("Pancreas", "PANC", names(hcp_PP_df))
names(hcp_PP_df) = gsub("Prostate", "PROS", names(hcp_PP_df))
names(hcp_PP_df) = gsub("SmallIntestine", "SMA", names(hcp_PP_df))
# Add empty riskmod and interAge columns
hcp_PP_df$interAge = hcp_PP_df$riskmod = list(character(0))
# PanelPRO with modifiers
hcp_PP_mod_df = data.frame(hcp_PP_df, hcp_markers_PP)


# Drop individuals with missing values
hcp_BM_df = hcp_BM_df[!drop_idx,]
hcp_BM_mod_df = hcp_BM_mod_df[!drop_idx,]
hcp_PP_df = hcp_PP_df[!drop_idx,]
hcp_PP_mod_df = hcp_PP_mod_df[!drop_idx,]


# Separate giant data frames into lists of families
# Only include families that have a proband

# BayesMendel without modifiers
hcp_families_BM = lapply(unique(hcp_BM_df$FamID), function(famID){
  if (sum(hcp_BM_df$isProband[hcp_BM_df$FamID==famID])==1) {
    return(hcp_BM_df[hcp_BM_df$FamID==famID, -1])
  } else {
    return(NULL)
  }
})
names(hcp_families_BM) = unique(hcp_BM_df$FamID)
hcp_families_BM[sapply(hcp_families_BM, is.null)] = NULL

# BayesMendel with modifiers
hcp_families_BM_mod = lapply(unique(hcp_BM_mod_df$FamID), function(famID){
  if (sum(hcp_BM_mod_df$isProband[hcp_BM_mod_df$FamID==famID])==1) {
    return(hcp_BM_mod_df[hcp_BM_mod_df$FamID==famID, -1])
  } else {
    return(NULL)
  }
})
names(hcp_families_BM_mod) = unique(hcp_BM_mod_df$FamID)
hcp_families_BM_mod[sapply(hcp_families_BM_mod, is.null)] = NULL

# PanelPRO without modifiers
hcp_families_PP = lapply(unique(hcp_PP_df$FamID), function(famID){
  if (sum(hcp_PP_df$isProband[hcp_PP_df$FamID==famID])==1) {
    return(hcp_PP_df[hcp_PP_df$FamID==famID, -1])
  } else {
    return(NULL)
  }
})
names(hcp_families_PP) = unique(hcp_PP_df$FamID)
hcp_families_PP[sapply(hcp_families_PP, is.null)] = NULL

# PanelPRO with modifiers
hcp_families_PP_mod = lapply(unique(hcp_PP_mod_df$FamID), function(famID){
  if (sum(hcp_PP_mod_df$isProband[hcp_PP_mod_df$FamID==famID])==1) {
    return(hcp_PP_mod_df[hcp_PP_mod_df$FamID==famID, -1])
  } else {
    return(NULL)
  }
})
names(hcp_families_PP_mod) = unique(hcp_PP_mod_df$FamID)
hcp_families_PP_mod[sapply(hcp_families_PP_mod, is.null)] = NULL


# Pull out mutation carrier status for each family
hcp_mutations = lapply(hcp_families_PP, function(fam){
  fam[mutations]
})

# Hack: because you dropped individuals with missing info, some individuals 
# will have parents with IDs that are no longer in the family, so replace 
# them with zeros
# Do appropriate processing for marker and germline testing

# BayesMendel without modifiers
hcp_families_BM = lapply(hcp_families_BM, function(fam){
  fam$MotherID[!sapply(fam$MotherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  fam$FatherID[!sapply(fam$FatherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  return(fam)
})

# BayesMendel with modifiers
hcp_families_BM_mod = lapply(hcp_families_BM_mod, function(fam){
  fam$MotherID[!sapply(fam$MotherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  fam$FatherID[!sapply(fam$FatherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  
  brca.marker.testing = fam[c("ER", "CK14", "CK5.6", "PR", "HER2")]
  if (all(rowSums(brca.marker.testing)==0)) {
    brca.marker.testing = NULL
  }
  mmr.marker.testing = data.frame(MSI=fam$MSI, location=0)
  if (all(rowSums(mmr.marker.testing)==0)) {
    mmr.marker.testing = NULL
  }
  
  brca.germline.testing = data.frame(fam[c("BRCA1", "BRCA2")], TestOrder=0)
  mmr.germline.testing = data.frame(fam[c("MLH1", "MSH2", "MSH6")], TestOrder=0)
  mela.germline.testing = data.frame(P16=fam$CDKN2A, TestOrder=0)
  proband_date = fam$germline.date[fam$isProband==1]
  if (!is.na(proband_date)) {
    brca.germline.testing[is.na(fam$germline.date) | 
                            fam$germline.date >= proband_date, 
                          c("BRCA1", "BRCA2")] = 0
    if (all(rowSums(brca.germline.testing)==0)) {
      brca.germline.testing = NULL
    }
    mmr.germline.testing[is.na(fam$germline.date) | 
                           fam$germline.date >= proband_date, 
                         c("MLH1", "MSH2", "MSH6")] = 0
    if (all(rowSums(mmr.germline.testing)==0)) {
      mmr.germline.testing = NULL
    }
    mela.germline.testing[is.na(fam$germline.date) | 
                            fam$germline.date >= proband_date, 
                          "P16"] = 0
    if (all(rowSums(mela.germline.testing)==0)) {
      mela.germline.testing = NULL
    }
  } else {
    brca.germline.testing = NULL
  }
  
  return(list(fam=fam, brcapro=list(germline.testing=brca.germline.testing, 
                                    marker.testing=brca.marker.testing), 
              mmrpro=list(germline.testing=mmr.germline.testing, 
                           marker.testing=mmr.marker.testing), 
              melapro=list(germline.testing=mela.germline.testing)))
})

# PanelPRO without modifiers
hcp_families_PP = lapply(hcp_families_PP, function(fam){
  fam$MotherID[!sapply(fam$MotherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  fam$FatherID[!sapply(fam$FatherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  fam$germline.date = NULL
  return(fam)
})

# PanelPRO with modifiers
hcp_families_PP_mod = lapply(hcp_families_PP_mod, function(fam){
  fam$MotherID[!sapply(fam$MotherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  fam$FatherID[!sapply(fam$FatherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  
  proband_date = fam$germline.date[fam$isProband==1]
  if (!is.na(proband_date)) {
    fam[is.na(fam$germline.date) | fam$germline.date >= proband_date, 
        mutations] = NA
  } else {
    fam[,mutations] = NA
  }
  fam$germline.date = NULL
  
  return(fam)
})


# 1218 families left
# Families without modifiers
save(hcp_families_BM, hcp_families_PP, 
     hcp_mutations, 
     file="hcp_families.rData")
# Families with modifiers
save(hcp_families_BM_mod, hcp_families_PP_mod, 
     hcp_mutations, 
     file="hcp_families_mod.rData")
