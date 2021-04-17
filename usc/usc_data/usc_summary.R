# PanelPRO-5BC after pre-processing
load("pp5/usc_families_mod.rData")
# Average number of people per family
mean(sapply(usc_families_PP, nrow))

# Number of probands who are carriers
colSums(sapply(c("ATM", "BRCA1", "BRCA2", 
                 "CHEK2", "PALB2"), function(mut) {
  sapply(usc_families_PP, function(fam){
    fam[[mut]][fam$isProband==1]
  })
}))
# Total number of probands who are carriers
sum(sapply(usc_families_PP, function(fam){
  any(fam[, c("ATM", "BRCA1", "BRCA2", 
              "CHEK2", "PALB2")][fam$isProband==1,] == 1)
}))
# Number of probands who are noncarriers
sum(sapply(usc_families_PP, function(fam){
  all(fam[, c("ATM", "BRCA1", "BRCA2", 
              "CHEK2", "PALB2")][fam$isProband==1,] == 0)
}))

# Number of probands affected for cancers in the model
colSums(sapply(c("BC", "OC"), function(canc) {
  sapply(usc_families_PP, function(fam){
    fam[[paste0("isAff", canc)]][fam$isProband==1]
  })
}))
# Number of probands affected for other cancers
sum(sapply(usc_families_PP, function(fam){
  fam[["OtherCancers"]][fam$isProband==1]
}), na.rm = TRUE)
# Total number of probands who are affected for cancer
sum(apply(cbind(
  sapply(c("BC", "OC"), function(canc) {
    sapply(usc_families_PP, function(fam){
      fam[[paste0("isAff", canc)]][fam$isProband==1]
    })
  }), 
  sapply(usc_families_PP, function(fam){
    fam[["OtherCancers"]][fam$isProband==1]
  })) == 1, 1, any))

# Number of carriers across all relatives
colSums(sapply(c("ATM", "BRCA1", "BRCA2", 
                 "CHEK2", "PALB2"), function(mut) {
                   sapply(usc_families_PP, function(fam){
                     sum(fam[[mut]])
                   })
                 }))
# Number of cancer cases across all relatives
colSums(sapply(c("BC", "OC"), function(canc) {
  sapply(usc_families_PP, function(fam){
    sum(fam[[paste0("isAff", canc)]])
  })
}))


# Number of families with at least one set of tumor marker results
colSums(sapply(c("ER", "CK14", "CK5.6", "PR", "HER2"), function(marker) {
  sapply(usc_families_PP_mod, function(fam){
    sum(fam[[marker]], na.rm=TRUE) > 0
  })
}))



# PanelPRO-11 after pre-processing
load("pp11/usc_families_mod.rData")
# Average number of people per family
mean(sapply(usc_families_PP, nrow))

# Number of probands who are carriers
colSums(sapply(c("ATM", "BRCA1", "BRCA2", "CDKN2A", 
                 "CHEK2", "EPCAM", "MLH1", "MSH2", 
                 "MSH6", "PALB2", "PMS2"), function(mut) {
                   sapply(usc_families_PP, function(fam){
                     fam[[mut]][fam$isProband==1]
                   })
                 }))
# Total number of probands who are carriers
sum(sapply(usc_families_PP, function(fam){
  any(fam[, c("ATM", "BRCA1", "BRCA2", "CDKN2A", 
              "CHEK2", "EPCAM", "MLH1", "MSH2", 
              "MSH6", "PALB2", "PMS2")][fam$isProband==1,] == 1)
}))
# Number of probands who are noncarriers
sum(sapply(usc_families_PP, function(fam){
  all(fam[, c("ATM", "BRCA1", "BRCA2", "CDKN2A", 
              "CHEK2", "EPCAM", "MLH1", "MSH2", 
              "MSH6", "PALB2", "PMS2")][fam$isProband==1,] == 0)
}))

# Number of probands affected for cancers in the model
colSums(sapply(c("BRA", "BC", "COL", "ENDO", 
                 "GAS", "KID", "MELA", "OC", 
                 "PANC", "PROS", "SMA"), function(canc) {
  sapply(usc_families_PP, function(fam){
    fam[[paste0("isAff", canc)]][fam$isProband==1]
  })
}))
# Number of probands affected for other cancers
sum(sapply(usc_families_PP, function(fam){
  fam[["OtherCancers"]][fam$isProband==1]
}), na.rm = TRUE)
# Total number of probands who are affected for cancer
sum(apply(cbind(
  sapply(c("BRA", "BC", "COL", "ENDO", 
           "GAS", "KID", "MELA", "OC", 
           "PANC", "PROS", "SMA"), function(canc) {
    sapply(usc_families_PP, function(fam){
      fam[[paste0("isAff", canc)]][fam$isProband==1]
    })
  }), 
  sapply(usc_families_PP, function(fam){
    fam[["OtherCancers"]][fam$isProband==1]
  })) == 1, 1, any))

# Number of carriers across all relatives
colSums(sapply(c("ATM", "BRCA1", "BRCA2", "CDKN2A", 
                 "CHEK2", "EPCAM", "MLH1", "MSH2", 
                 "MSH6", "PALB2", "PMS2"), function(mut) {
                   sapply(usc_families_PP, function(fam){
                     sum(fam[[mut]])
                   })
                 }))
# Number of cancer cases across all relatives
colSums(sapply(c("BRA", "BC", "COL", "ENDO", 
                 "GAS", "KID", "MELA", "OC", 
                 "PANC", "PROS", "SMA"), function(canc) {
  sapply(usc_families_PP, function(fam){
    sum(fam[[paste0("isAff", canc)]])
  })
}))

# Allele frequencies among probands
colSums(sapply(c("ATM", "BRCA1", "BRCA2", "CDKN2A", 
                 "CHEK2", "EPCAM", "MLH1", "MSH2", 
                 "MSH6", "PALB2", "PMS2"), function(mut) {
                   sapply(usc_families_PP, function(fam){
                     fam[[mut]][fam$isProband==1]
                   })
                 })) / length(usc_families_PP) / 2
# Allele frequencies among all relatives (not all members were tested, though)
colSums(sapply(c("ATM", "BRCA1", "BRCA2", "CDKN2A", 
                 "CHEK2", "EPCAM", "MLH1", "MSH2", 
                 "MSH6", "PALB2", "PMS2"), function(mut) {
                   sapply(usc_families_PP, function(fam){
                     sum(fam[[mut]])
                   })
                 })) / 
  sum(sapply(usc_families_PP, nrow)) / 2


# Number of families with at least one set of tumor marker results
colSums(sapply(c("ER", "CK14", "CK5.6", "PR", "HER2", "MSI"), function(marker) {
  sapply(usc_families_PP_mod, function(fam){
    sum(fam[[marker]], na.rm=TRUE) > 0
  })
}))
