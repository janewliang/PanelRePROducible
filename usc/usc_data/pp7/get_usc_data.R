# Load USC data
load("../usc_dat.RData")

# 6 mutations to consider (PANC not included)
mutations = c("BRCA1", "BRCA2", "CDKN2A", "MLH1", "MSH2", "MSH6")

# Drop probands who are carriers of other mutations, but not one of the 6
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
# Noncarriers = negative for all mutations in model
# Drop probands with VUSs and no pathogenic genes in model
isVUS = 
  # Individual has a VUS
  apply((usc[,paste0("res", 1:5)] == 0) & # VUS
          apply(usc[,paste0("gene", 1:5)], 1:2, function(x){ # One of the genes
            x %in% mutations}), 
        1, any) & 
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

# 6 cancers to consider
# Create data frame of cancers
usc_cancers = usc[,c("AffectedBreast", "AgeBreast", 
                     "AffectedColon", "AgeColon", 
                     "AffectedEndometrium", "AgeEndometrium", 
                     "AffectedSkin", "AgeSkin", 
                     "AffectedOvary", "AgeOvary", 
                     "AffectedPancreas", "AgePancreas")]


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

# Set maximum age to PanelPRO's MAXAGE (94)
usc_df$AgeDeath[usc_df$AgeDeath > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeBreast[usc_df$AgeBreast > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeColon[usc_df$AgeColon > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeEndometrium[usc_df$AgeEndometrium > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeSkin[usc_df$AgeSkin > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgeOvary[usc_df$AgeOvary > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
usc_df$AgePancreas[usc_df$AgePancreas > PanelPRO:::MAXAGE] = PanelPRO:::MAXAGE
# Hack to substitute affected age=1 with 2
# BM models use 1 to code for unknown affected age
usc_df$AgeBreast[usc_df$AgeBreast==1 & usc_df$AffectedBreast==1] = 2
usc_df$AgeOvary[usc_df$AgeOvary==1 & usc_df$AffectedOvary==1] = 2
usc_df$AgeColon[usc_df$AgeColon==1 & usc_df$AffectedColon==1] = 2
usc_df$AgeEndometrium[usc_df$AgeEndometrium==1 & usc_df$AffectedEndometrium==1] = 2
usc_df$AgeSkin[usc_df$AgeSkin==1 & usc_df$AffectedSkin==1] = 2


# BayesMendel
usc_BM_df = data.frame(usc_df, AgeBreastContralateral=0)

# PanelPRO
usc_PP_df = data.frame(usc_df, race="All_Races", stringsAsFactors=FALSE)
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


# Drop individuals with missing values
usc_BM_df = usc_BM_df[!drop_idx,]
usc_PP_df = usc_PP_df[!drop_idx,]


# Separate giant data frames into lists of families
# Only include families that have a proband

# BayesMendel
usc_families_BM = lapply(unique(usc_BM_df$FamID), function(famID){
  if (sum(usc_BM_df$isProband[usc_BM_df$FamID==famID])==1) {
    return(usc_BM_df[usc_BM_df$FamID==famID, -1])
  } else {
    return(NULL)
  }
})
names(usc_families_BM) = unique(usc_BM_df$FamID)
usc_families_BM[sapply(usc_families_BM, is.null)] = NULL

# PanelPRO
usc_families_PP = lapply(unique(usc_PP_df$FamID), function(famID){
  if (sum(usc_PP_df$isProband[usc_PP_df$FamID==famID])==1) {
    return(usc_PP_df[usc_PP_df$FamID==famID, -1])
  } else {
    return(NULL)
  }
})
names(usc_families_PP) = unique(usc_PP_df$FamID)
usc_families_PP[sapply(usc_families_PP, is.null)] = NULL


# Pull out mutation carrier status for each family
usc_mutations = lapply(usc_families_PP, function(fam){
  fam[mutations]
})

# Hack: because you dropped individuals with missing info, some individuals 
# will have parents with IDs that are no longer in the family, so replace 
# them with zeros
# Do appropriate processing for marker and germline testing

# BayesMendel
usc_families_BM = lapply(usc_families_BM, function(fam){
  fam$MotherID[!sapply(fam$MotherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  fam$FatherID[!sapply(fam$FatherID, function(x){
    x %in% c(0, fam$ID)
  })] = 0
  return(fam)
})

# PanelPRO
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


# 1688 families left
save(usc_families_BM, usc_families_PP, 
     usc_mutations, 
     file="usc_families.rData")
