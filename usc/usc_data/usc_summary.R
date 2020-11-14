# Load USC data
load("usc_dat_pp11.RData")

# Subset just the probands
proband_usc = usc[usc$Proband==1,]

# The mutation carriers
proband_gene = proband_usc[,paste0("gene", 1:5)]
# Mutation status (carrier, VUS, etc.)
proband_res = proband_usc[,paste0("res", 1:5)]

# Get a table of the mutation-positive carriers (probands only)
table(proband_gene[proband_res==1])

# Table of proband affected status
colSums(proband_usc[,c("AffectedBrain", 
                       "AffectedBreast", 
                       "AffectedColon", 
                       "AffectedEndometrium", 
                       "AffectedGastric", 
                       "AffectedKidney", 
                       "AffectedSkin", 
                       "AffectedOvary", 
                       "AffectedPancreas", 
                       "AffectedProstate", 
                       "AffectedSmallIntestine")])



# PanelPRO-5BC after pre-processing
load("pp5/usc_families.rData")
colSums(sapply(c("ATM", "BRCA1", "BRCA2", "CHEK2", "PALB2"), function(mut) {
  sapply(usc_families_PP, function(fam){
    fam[[mut]][fam$isProband==1]
  })
}))
colSums(sapply(c("BC", "OC"), function(canc) {
  sapply(usc_families_PP, function(fam){
    fam[[paste0("isAff", canc)]][fam$isProband==1]
  })
}))


# PanelPRO-11 after pre-processing
load("pp11/usc_families.rData")
colSums(sapply(c("ATM", "BRCA1", "BRCA2", "CDKN2A", 
                 "CHEK2", "EPCAM", "MLH1", "MSH2", 
                 "MSH6", "PALB2", "PMS2"), function(mut) {
  sapply(usc_families_PP, function(fam){
    fam[[mut]][fam$isProband==1]
  })
}))
colSums(sapply(c("BRA", "BC", "COL", "ENDO", "GAS", "KID", 
                 "MELA", "OC", "PANC", "PROS", "SMA"), function(canc) {
  sapply(usc_families_PP, function(fam){
    fam[[paste0("isAff", canc)]][fam$isProband==1]
  })
}))


# PanelPRO-5BC risk modifiers after pre-processing
load("pp5/usc_families_mod.rData")
colSums(sapply(c("ER", "CK14", "CK5.6", "PR", "HER2"), function(marker) {
  sapply(usc_families_PP_mod, function(fam){
    sum(fam[[marker]], na.rm=TRUE) > 0
  })
}))


# PanelPRO-11 risk modifiers after pre-processing
load("pp11/usc_families_mod.rData")
colSums(sapply(c("ER", "CK14", "CK5.6", "PR", "HER2", "MSI"), function(marker) {
  sapply(usc_families_PP_mod, function(fam){
    sum(fam[[marker]], na.rm=TRUE) > 0
  })
}))
