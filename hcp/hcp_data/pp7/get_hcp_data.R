# Load HCP data
load("../hcp_data.RData")

# 6 mutations to consider (PANC not included)
mutations = c("BRCA1", "BRCA2", "CDKN2A", "MLH1", "MSH2", "MSH6")

# Drop probands who are carriers of other mutations, but not one of the 6
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

# 6 cancers to consider
# Create data frame of cancers
hcp_cancers = hcp[,c("AffectedBreast", "AgeBreast", 
                     "AffectedColon", "AgeColon", 
                     "AffectedEndometrium", "AgeEndometrium", 
                     "AffectedSkin", "AgeSkin", 
                     "AffectedOvary", "AgeOvary", 
                     "AffectedPancreas", "AgePancreas")]


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

# Set maximum age to PanelPRO's MAXAGE (94)
hcp_df$AgeDeath[hcp_df$AgeDeath > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeBreast[hcp_df$AgeBreast > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeColon[hcp_df$AgeColon > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeEndometrium[hcp_df$AgeEndometrium > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeSkin[hcp_df$AgeSkin > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgeOvary[hcp_df$AgeOvary > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
hcp_df$AgePancreas[hcp_df$AgePancreas > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
# Hack to substitute affected age=1 with 2
# BM models use 1 to code for unknown affected age
hcp_df$AgeBreast[hcp_df$AgeBreast==1 & hcp_df$AffectedBreast==1] = 2
hcp_df$AgeOvary[hcp_df$AgeOvary==1 & hcp_df$AffectedOvary==1] = 2
hcp_df$AgeColon[hcp_df$AgeColon==1 & hcp_df$AffectedColon==1] = 2
hcp_df$AgeEndometrium[hcp_df$AgeEndometrium==1 & hcp_df$AffectedEndometrium==1] = 2
hcp_df$AgeSkin[hcp_df$AgeSkin==1 & hcp_df$AffectedSkin==1] = 2


# BayesMendel
hcp_BM_df = data.frame(hcp_df, AgeBreastContralateral=0)

# PanelPRO
hcp_PP_df = data.frame(hcp_df, race="All_Races", stringsAsFactors=FALSE)
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


# Drop individuals with missing values
hcp_BM_df = hcp_BM_df[!drop_idx,]
hcp_PP_df = hcp_PP_df[!drop_idx,]


# Separate giant data frames into lists of families
# Only include families that have a proband

# BayesMendel
hcp_families_BM = lapply(unique(hcp_BM_df$FamID), function(famID){
  if (sum(hcp_BM_df$isProband[hcp_BM_df$FamID==famID])==1) {
    return(hcp_BM_df[hcp_BM_df$FamID==famID, -1])
  } else {
    return(NULL)
  }
})
names(hcp_families_BM) = unique(hcp_BM_df$FamID)
hcp_families_BM[sapply(hcp_families_BM, is.null)] = NULL

# PanelPRO
hcp_families_PP = lapply(unique(hcp_PP_df$FamID), function(famID){
  if (sum(hcp_PP_df$isProband[hcp_PP_df$FamID==famID])==1) {
    return(hcp_PP_df[hcp_PP_df$FamID==famID, -1])
  } else {
    return(NULL)
  }
})
names(hcp_families_PP) = unique(hcp_PP_df$FamID)
hcp_families_PP[sapply(hcp_families_PP, is.null)] = NULL


# Pull out mutation carrier status for each family
hcp_mutations = lapply(hcp_families_PP, function(fam){
  fam[mutations]
})

# Hack: because you dropped individuals with missing info, some individuals 
# will have parents with IDs that are no longer in the family, so replace 
# them with zeros
# Do appropriate processing for marker and germline testing

# BayesMendel
hcp_families_BM = lapply(hcp_families_BM, function(fam){
  fam$MotherID[!sapply(fam$MotherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  fam$FatherID[!sapply(fam$FatherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  return(fam)
})

# PanelPRO
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


# 1169 families left
save(hcp_families_BM, hcp_families_PP, 
     hcp_mutations, 
     file="hcp_families.rData")
