# PanelPRO-5BC families after pre-processing
load("pp5/hcp_families.rData")
colSums(sapply(c("ATM", "BRCA1", "BRCA2", "CHEK2", "PALB2"), function(mut) {
  sapply(hcp_families_PP, function(fam){
    fam[[mut]][fam$isProband==1]
  })
}))
colSums(sapply(c("BC", "OC"), function(canc) {
  sapply(hcp_families_PP, function(fam){
    fam[[paste0("isAff", canc)]][fam$isProband==1]
  })
}))


# PanelPRO-11 families after pre-processing
load("pp11/hcp_families.rData")
colSums(sapply(c("ATM", "BRCA1", "BRCA2", "CDKN2A", 
                 "CHEK2", "EPCAM", "MLH1", "MSH2", 
                 "MSH6", "PALB2", "PMS2"), function(mut) {
  sapply(hcp_families_PP, function(fam){
    fam[[mut]][fam$isProband==1]
  })
}))
colSums(sapply(c("BRA", "BC", "COL", "ENDO", "GAS", "KID", 
                 "MELA", "OC", "PANC", "PROS", "SMA"), function(canc) {
  sapply(hcp_families_PP, function(fam){
    fam[[paste0("isAff", canc)]][fam$isProband==1]
  })
}))


# PanelPRO-5BC risk modifiers after pre-processing
load("pp5/hcp_families_mod.rData")
colSums(sapply(c("ER", "CK14", "CK5.6", "PR", "HER2"), function(marker) {
  sapply(hcp_families_PP_mod, function(fam){
    sum(fam[[marker]], na.rm=TRUE) > 0
  })
}))


# PanelPRO-11 risk modifiers after pre-processing
load("pp11/hcp_families_mod.rData")
colSums(sapply(c("ER", "CK14", "CK5.6", "PR", "HER2", "MSI"), function(marker) {
  sapply(hcp_families_PP_mod, function(fam){
    sum(fam[[marker]], na.rm=TRUE) > 0
  })
}))
