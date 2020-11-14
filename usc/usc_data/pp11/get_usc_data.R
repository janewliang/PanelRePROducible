# Load USC data
load("../usc_dat_pp11.RData")

# 11 mutations to consider
mutations = c("ATM", "BRCA1", "BRCA2", "CDKN2A", "CHEK2", "EPCAM",
              "MLH1", "MSH2", "MSH6", "PALB2", "PMS2")

# Drop probands who are carriers of other mutations, but not one of the 11
isCarrOtherMut = apply((usc[,paste0("res", 1:5)] == 1) & 
                         !(apply(usc[,paste0("gene", 1:5)], 1:2, function(x){
                           x %in% c(mutations, NA)})), 
                       1, any) &
  !(apply((usc[,paste0("res", 1:5)] == 1) & 
            apply(usc[,paste0("gene", 1:5)], 1:2, function(x){
              x %in% mutations}), 
          1, any))
isCarrOtherMut[is.na(isCarrOtherMut)==TRUE] = FALSE
usc = usc[!(usc$Proband==1 & isCarrOtherMut==TRUE),]

# Carriers = at least one pathogenic mutation
# Noncarriers = negative for all mutations
# Drop probands with VUSs and no pathogenic genes
isVUS = 
  # Individual has VUSs
  apply(usc[,paste0("res", 1:5)], 1, function(x){ 0 %in% x}) & 
  # Individual has no pathogenic genes
  !apply((usc[,paste0("res", 1:5)] == 1) & # Pathogenic
          apply(usc[,paste0("gene", 1:5)], 1:2, function(x){ # One of the genes
            x %in% mutations}), 
        1, any)
usc = usc[!(usc$Proband==1 & isVUS==TRUE),]

# Create data frame of mutations
usc_muts = data.frame(sapply(mutations, function(mut){
  (usc$gene1==mut & usc$res1==1) | (usc$gene2==mut & usc$res2==1) | 
    (usc$gene3==mut & usc$res3==1) | (usc$gene4==mut & usc$res4==1) | 
    (usc$gene5==mut & usc$res5==1)
}))
usc_muts[is.na(usc_muts)] = 0

# 11 cancers to consider
# Create data frame of cancers
usc_cancers = usc[,c("AffectedBrain", "AgeBrain", 
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
usc_markers_BM = usc[markers]
# Tumor markers for PanelPRO 
usc_markers_PP = usc_markers_BM
# Set 0 to NA (not tested)
usc_markers_PP[usc_markers_PP==0] = NA
# Set 2 to 0 (negative test)
usc_markers_PP[usc_markers_PP==2] = 0


# Build data frame of all families
usc_df = data.frame(FamID=usc$FamID, 
                    ID = usc$ID, 
                    MotherID = usc$MotherID, 
                    FatherID = usc$FatherID, 
                    Gender = usc$Gender, 
                    isProband = usc$Proband, 
                    Twins = 0, 
                    ethnic = "nonAJ",
                    AgeDeath = as.numeric(usc$Age.of.entry), 
                    usc_cancers, 
                    Death = as.numeric(usc$Deceased.Status), 
                    usc_muts, stringsAsFactors=FALSE)
# Identify individuals with missing values
drop_idx = rowSums(is.na(usc_df)) > 0
# Date for germline testing
germline.date = as.Date(usc$GenTestTable.Result.Date, format="%Y-%m-%d")

# Set maximum age to PanelPRO's MAXAGE (94)
usc_df$AgeDeath[usc_df$AgeDeath > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeBrain[usc_df$AgeBrain > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeBreast[usc_df$AgeBreast > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeColon[usc_df$AgeColon > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeEndometrium[usc_df$AgeEndometrium > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeGastric[usc_df$AgeGastric > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeKidney[usc_df$AgeKidney > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeSkin[usc_df$AgeSkin > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeOvary[usc_df$AgeOvary > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgePancreas[usc_df$AgePancreas > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeProstate[usc_df$AgeProstate > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeSmallIntestine[usc_df$AgeSmallIntestine > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
# Hack to substitute affected age=1 with 2
# BM models use 1 to code for unknown affected age
usc_df$AgeBreast[usc_df$AgeBreast==1 & usc_df$AffectedBreast==1] = 2
usc_df$AgeOvary[usc_df$AgeOvary==1 & usc_df$AffectedOvary==1] = 2
usc_df$AgeColon[usc_df$AgeColon==1 & usc_df$AffectedColon==1] = 2
usc_df$AgeEndometrium[usc_df$AgeEndometrium==1 & usc_df$AffectedEndometrium==1] = 2
usc_df$AgeSkin[usc_df$AgeSkin==1 & usc_df$AffectedSkin==1] = 2


# BayesMendel without modifiers
usc_BM_df = data.frame(usc_df, AgeBreastContralateral=0, 
                       germline.date=germline.date)
# BayesMendel with modifiers
usc_BM_mod_df = data.frame(usc_BM_df, usc_markers_BM)

# PanelPRO
usc_PP_df = data.frame(usc_df, race="All_Races", 
                       germline.date=germline.date, 
                       stringsAsFactors=FALSE)
# Re-name columns
names(usc_PP_df) = gsub("Gender", "Sex", names(usc_PP_df))
names(usc_PP_df) = gsub("AgeDeath", "CurAge", names(usc_PP_df))
names(usc_PP_df) = gsub("Death", "isDead", names(usc_PP_df))
names(usc_PP_df) = gsub("Affected", "isAff", names(usc_PP_df))
names(usc_PP_df) = gsub("ethnic", "Ancestry", names(usc_PP_df))
# Use cancer codes
names(usc_PP_df) = gsub("Brain", "BRA", names(usc_PP_df))
names(usc_PP_df) = gsub("Breast", "BC", names(usc_PP_df))
names(usc_PP_df) = gsub("Colon", "COL", names(usc_PP_df))
names(usc_PP_df) = gsub("Endometrium", "ENDO", names(usc_PP_df))
names(usc_PP_df) = gsub("Gastric", "GAS", names(usc_PP_df))
names(usc_PP_df) = gsub("Kidney", "KID", names(usc_PP_df))
names(usc_PP_df) = gsub("Skin", "MELA", names(usc_PP_df))
names(usc_PP_df) = gsub("Ovary", "OC", names(usc_PP_df))
names(usc_PP_df) = gsub("Pancreas", "PANC", names(usc_PP_df))
names(usc_PP_df) = gsub("Prostate", "PROS", names(usc_PP_df))
names(usc_PP_df) = gsub("SmallIntestine", "SMA", names(usc_PP_df))
# Add empty riskmod and interAge columns
usc_PP_df$interAge = usc_PP_df$riskmod = list(character(0))
# PanelPRO with modifiers
usc_PP_mod_df = data.frame(usc_PP_df, usc_markers_PP)


# Drop individuals with missing values
usc_BM_df = usc_BM_df[!drop_idx,]
usc_BM_mod_df = usc_BM_mod_df[!drop_idx,]
usc_PP_df = usc_PP_df[!drop_idx,]
usc_PP_mod_df = usc_PP_mod_df[!drop_idx,]


# Separate giant data frames into lists of families
# Only include families that have a proband

# BayesMendel without modifiers
usc_families_BM = lapply(unique(usc_BM_df$FamID), function(famID){
  if (sum(usc_BM_df$isProband[usc_BM_df$FamID==famID])==1) {
    return(usc_BM_df[usc_BM_df$FamID==famID, -1])
  } else {
    return(NULL)
  }
})
names(usc_families_BM) = unique(usc_BM_df$FamID)
usc_families_BM[sapply(usc_families_BM, is.null)] = NULL

# BayesMendel with modifiers
usc_families_BM_mod = lapply(unique(usc_BM_mod_df$FamID), function(famID){
  if (sum(usc_BM_mod_df$isProband[usc_BM_mod_df$FamID==famID])==1) {
    return(usc_BM_mod_df[usc_BM_mod_df$FamID==famID, -1])
  } else {
    return(NULL)
  }
})
names(usc_families_BM_mod) = unique(usc_BM_mod_df$FamID)
usc_families_BM_mod[sapply(usc_families_BM_mod, is.null)] = NULL

# PanelPRO without modifiers
usc_families_PP = lapply(unique(usc_PP_df$FamID), function(famID){
  if (sum(usc_PP_df$isProband[usc_PP_df$FamID==famID])==1) {
    return(usc_PP_df[usc_PP_df$FamID==famID, -1])
  } else {
    return(NULL)
  }
})
names(usc_families_PP) = unique(usc_PP_df$FamID)
usc_families_PP[sapply(usc_families_PP, is.null)] = NULL

# PanelPRO with modifiers
usc_families_PP_mod = lapply(unique(usc_PP_mod_df$FamID), function(famID){
  if (sum(usc_PP_mod_df$isProband[usc_PP_mod_df$FamID==famID])==1) {
    return(usc_PP_mod_df[usc_PP_mod_df$FamID==famID, -1])
  } else {
    return(NULL)
  }
})
names(usc_families_PP_mod) = unique(usc_PP_mod_df$FamID)
usc_families_PP_mod[sapply(usc_families_PP_mod, is.null)] = NULL


# Pull out mutation carrier status for each family
usc_mutations = lapply(usc_families_PP, function(fam){
  fam[mutations]
})

# Hack: because you dropped individuals with missing info, some individuals 
# will have parents with IDs that are no longer in the family, so replace 
# them with zeros
# Do appropriate processing for marker and germline testing

# BayesMendel without modifiers
usc_families_BM = lapply(usc_families_BM, function(fam){
  fam$MotherID[!sapply(fam$MotherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  fam$FatherID[!sapply(fam$FatherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  return(fam)
})

# BayesMendel with modifiers
usc_families_BM_mod = lapply(usc_families_BM_mod, function(fam){
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
usc_families_PP = lapply(usc_families_PP, function(fam){
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
usc_families_PP_mod = lapply(usc_families_PP_mod, function(fam){
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
save(usc_families_BM, usc_families_PP, 
     usc_mutations, 
     file="usc_families.rData")
# Families with modifiers
save(usc_families_BM_mod, usc_families_PP_mod, 
     usc_mutations, 
     file="usc_families_mod.rData")

# Family sizes
fam_size = sapply(usc_families_PP_mod, nrow)
table(fam_size) # 48 have only 1 person (proband)
