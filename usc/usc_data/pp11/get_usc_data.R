# Load USC data
load("../../usc_dat.RData")
usc_PP_df = usc

# 11 mutations to consider
mutations = c("ATM", "BRCA1", "BRCA2", "CDKN2A", "CHEK2", "EPCAM",
              "MLH1", "MSH2", "MSH6", "PALB2", "PMS2")
# Columns corresponding to germline testing results
gene_cols = 58:84
# Indicator for which germline testing columns correspond to model mutations
isMut = names(usc_PP_df)[gene_cols] %in% mutations

# Drop probands who are carriers of other mutations, but not one of the 24
isCarrOtherMut = 
  (apply(usc_PP_df[,gene_cols][,!isMut] == 1, 1, any)) & 
  !(apply(usc_PP_df[,gene_cols][,isMut] == 1, 1, any))
isCarrOtherMut[is.na(isCarrOtherMut)==TRUE] = FALSE
usc_PP_df = usc_PP_df[!(usc_PP_df$isProband==1 & isCarrOtherMut==TRUE),]

# Carriers = at least one pathogenic mutation
# Noncarriers = negative for all mutations in model
# Drop probands with VUSs and no pathogenic genes in model
isVUS = 
  (apply(usc_PP_df[,gene_cols][,isMut] == 2, 1, any)) & 
  !(apply(usc_PP_df[,gene_cols][,isMut] == 1, 1, any))
isVUS[is.na(isVUS)==TRUE] = FALSE
usc_PP_df = usc_PP_df[!(usc_PP_df$isProband==1 & isVUS==TRUE),]

# Add empty riskmod and interAge columns
usc_PP_df$interAge = usc_PP_df$riskmod = list(character(0))


# PanelPRO-formatted families
usc_families_PP = lapply(unique(usc_PP_df$FamID), function(famID){
  if (sum(usc_PP_df$isProband[usc_PP_df$FamID==famID])==1) {
    return(usc_PP_df[usc_PP_df$FamID==famID, ])
  } else {
    return(NULL)
  }
})
names(usc_families_PP) = unique(usc_PP_df$FamID)
usc_families_PP[sapply(usc_families_PP, is.null)] = NULL


# Extract germline testing results for probands and mutations of interest
proband_muts = t(sapply(usc_families_PP, function(fam){
  as.numeric(unlist(fam[fam$isProband==1,mutations]))
}))
colnames(proband_muts) = mutations

# Set VUSs to 0
proband_muts[proband_muts == 2] = 0 
# Convert to data frame
proband_muts = data.frame(proband_muts)

# 1506 families left
# Families without modifiers and with mutation status
save(usc_families_PP, proband_muts, file="usc_families.rData")
