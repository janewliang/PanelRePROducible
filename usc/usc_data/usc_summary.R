# Function to generate "Table 1" information
# usc_list = list of pedigrees in PanelPRO-ready format
# genes = genes in the model/to consider
# cancers = cancers in the model/to consider
# markers = tumor markers in the model/to consider (optional)
# Returns a list with up to 4 tables: 
# - fam_size_tab = summary statistics for the number of people in each family 
# - proband_carrier_tab = number of probands who are carriers
# - proband_cancer_tab = number of probands who are affected for cancer
# - marker_tab = number of families where tumor marker information is available 
get_table1 = function(usc_list, genes, cancers, markers = NULL) {
  
  # Summary statistics for the number of people in each family
  fam_size_tab = summary(sapply(usc_list, nrow))
  
  # Distribution of ancestry for probands
  proband_ancestry_tab = table(sapply(usc_list, function(fam) {
    fam$Ancestry[fam$isProband==1]
  }))
  
  # Distribution of race for probands
  proband_race_tab = table(sapply(usc_list, function(fam) {
    fam$race[fam$isProband==1]
  }))
  
  # Number of probands who are carriers
  proband_carrier_tab = 
    c(# Number of probands who are carriers for genes in the model
      colSums(sapply(genes, function(mut) {
        sapply(usc_list, function(fam){
          fam[[mut]][fam$isProband==1]
        })
      })), 
      # Total number of probands who are carriers
      AnyCarrier = sum(sapply(usc_list, function(fam){
        any(fam[, genes][fam$isProband==1,] == 1)
      })), 
      # Number of probands who are noncarriers
      NonCarrier = sum(sapply(usc_list, function(fam){
        all(fam[, genes][fam$isProband==1,] == 0)
      }))
      )
  
  # ALl cancer columns indices
  all_cancer_col_idx = c(grep("isAff*", names(usc_list[[1]])), 
                         which(names(usc_list[[1]]) == "OtherCancerCount"))
  # Column indices for other cancers
  other_cancer_col_idx = 
    all_cancer_col_idx[!(all_cancer_col_idx %in%
                           sapply(cancers, function(canc) {
                             which(names(usc_list[[1]]) == paste0("isAff", canc))
                           }))]
  
  # Number of probands who are affected for cancer
  proband_cancer_tab = 
    c(# Number of probands affected for cancers in the model
      colSums(sapply(cancers, function(canc) {
                         sapply(usc_list, function(fam){
                           fam[[paste0("isAff", canc)]][fam$isProband==1]
                         })
                       })), 
      # Number of probands affected for other cancers
      OtherCancer = sum(sapply(usc_list, function(fam){
        any(fam[fam$isProband==1,other_cancer_col_idx] > 0)
      }), na.rm = TRUE), 
      # Total number of probands who are affected for cancer
      AnyCancer = sum(sapply(usc_list, function(fam){
        any(fam[fam$isProband==1,all_cancer_col_idx] > 0)
      }), na.rm = TRUE), 
      # Number of probands who are unaffected for cancer
      NoCancer = sum(sapply(usc_list, function(fam){
        all(fam[fam$isProband==1,all_cancer_col_idx] == 0)
      }))
      )
  
  if (!is.null(markers)) {
    # Number of families where tumor marker information is available 
    marker_tab = c(
      colSums(sapply(markers, function(marker) {
        sapply(usc_list, function(fam){
          sum(fam[[marker]], na.rm=TRUE) > 0
        })
      })), 
      # Total number of families where tumor marker information is available 
      AnyMarker = sum(sapply(usc_list, function(fam){
        any(sum(fam[, markers], na.rm=TRUE) > 0)
      }), na.rm = TRUE)
    )
    
    # Return list of tables
    return(list(fam_size_tab = fam_size_tab, 
                proband_ancestry_tab = proband_ancestry_tab,
                proband_race_tab = proband_race_tab, 
                proband_carrier_tab = proband_carrier_tab, 
                proband_cancer_tab = proband_cancer_tab, 
                marker_tab = marker_tab))
  } else {
    # Return list of tables
    return(list(fam_size_tab = fam_size_tab, 
                proband_ancestry_tab = proband_ancestry_tab,
                proband_race_tab = proband_race_tab, 
                proband_carrier_tab = proband_carrier_tab, 
                proband_cancer_tab = proband_cancer_tab))
  }
  
}


# PanelPRO-5BC
load("pp5/usc_families.rData")
load("../../usc/pp5/results/diagnostics/diagnostics.rData")
usc_families_PP = usc_families_PP[!failed_families]
pp5_tables = get_table1(usc_families_PP, 
                        genes = c("ATM", "BRCA1", "BRCA2", 
                                  "CHEK2", "PALB2"), 
                        cancers = c("BC", "OC"), 
                        markers = c("ER", "CK14", "CK5.6", 
                                    "PR", "HER2"))


# PanelPRO-11
load("pp11/usc_families.rData")
load("../../usc/pp11/results/diagnostics/diagnostics.rData")
usc_families_PP = usc_families_PP[!failed_families]
pp11_tables = get_table1(usc_families_PP, 
                         genes = c("ATM", "BRCA1", "BRCA2", "CDKN2A", 
                                   "CHEK2", "EPCAM", "MLH1", "MSH2", 
                                   "MSH6", "PALB2", "PMS2"), 
                         cancers = c("BRA", "BC", "COL", "ENDO", 
                                     "GAS", "KID", "MELA", "OC", 
                                     "PANC", "PROS", "SI"), 
                         markers = c("ER", "CK14", "CK5.6", 
                                     "PR", "HER2", "MSI"))


# PanelPRO-24 after pre-processing
load("pp24/usc_families.rData")
load("../../usc/pp24/results/diagnostics/diagnostics.rData")
usc_families_PP = usc_families_PP[!failed_families]
pp24_tables = get_table1(usc_families_PP, 
                         genes = c("APC", "ATM", "BARD1", "BMPR1A", 
                                   "BRCA1", "BRCA2", "BRIP1", "CDH1", 
                                   "CDK4", "CDKN2A", "CHEK2", "EPCAM", 
                                   "MLH1", "MSH2", "MSH6", "MUTYH", 
                                   "NBN", "PALB2", "PMS2", "PTEN", 
                                   "RAD51C", "RAD51D", "STK11", "TP53"), 
                         cancers = c("BRA", "BC", "COL", "ENDO", 
                                     "GAS", "KID", "MELA", "OC", 
                                     "PANC", "PROS", "SI"), 
                         markers = c("ER", "CK14", "CK5.6", 
                                     "PR", "HER2", "MSI"))
