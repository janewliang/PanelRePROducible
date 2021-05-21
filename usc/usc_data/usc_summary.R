# Function to generate "Table 1" information
# usc_list = list of pedigrees in PanelPRO-ready format
# usc_list_mod = list of pedigrees in PanelPRO-ready format with risk- 
# modifying information
# genes = genes in the model/to consider
# cancers = cancers in the model/to consider
# markers = tumor markers in the model/to consider (optional)
# Returns a list with up to 4 tables: 
# - fam_size_tab = summary statistics for the number of people in each family 
# - proband_carrier_tab = number of probands who are carriers
# - proband_cancer_tab = number of probands who are affected for cancer
# - marker_tab = number of families where tumor marker information is available 
get_table1 = function(usc_list, usc_list_mod, 
                      genes, cancers, markers = NULL) {
  
  # Summary statistics for the number of people in each family
  fam_size_tab = summary(sapply(usc_list, nrow))
  
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
        fam[["OtherCancers"]][fam$isProband==1]
      }), na.rm = TRUE), 
      # Total number of probands who are affected for cancer
      AnyCancer = sum(apply(cbind(
        sapply(cancers, function(canc) {
          sapply(usc_list, function(fam){
            fam[[paste0("isAff", canc)]][fam$isProband==1]
          })
        }), 
        sapply(usc_families_PP, function(fam){
          fam[["OtherCancers"]][fam$isProband==1]
        })) == 1, 1, any)), 
      # Number of probands who are unaffected for cancer
      NoCancer = sum(sapply(usc_list, function(fam){
        all(fam[, c(paste0("isAff", cancers), "OtherCancers")][fam$isProband==1,] == 0)
      }))
    )
  
  if (!is.null(markers)) {
    # Number of families where tumor marker information is available 
    marker_tab = colSums(sapply(markers, function(marker) {
      sapply(usc_list_mod, function(fam){
        sum(fam[[marker]], na.rm=TRUE) > 0
      })
    }))
    
    # Return list of tables
    return(list(fam_size_tab = fam_size_tab, 
                proband_carrier_tab = proband_carrier_tab, 
                proband_cancer_tab = proband_cancer_tab, 
                marker_tab = marker_tab))
  } else {
    # Return list of tables
    return(list(fam_size_tab = fam_size_tab, 
                carrier_tab = carrier_tab, 
                cancer_tab = cancer_tab))
  }
  
}


# PanelPRO-5BC
load("pp5/usc_families.rData")
load("pp5/usc_families_mod.rData")
pp5_tables = get_table1(usc_families_PP, usc_families_PP_mod, 
                        genes = c("ATM", "BRCA1", "BRCA2", 
                                  "CHEK2", "PALB2"), 
                        cancers = c("BC", "OC"), 
                        markers = c("ER", "CK14", "CK5.6", 
                                    "PR", "HER2"))


# PanelPRO-11
load("pp11/usc_families.rData")
load("pp11/usc_families_mod.rData")
pp11_tables = get_table1(usc_families_PP, usc_families_PP_mod, 
                         genes = c("ATM", "BRCA1", "BRCA2", "CDKN2A", 
                                   "CHEK2", "EPCAM", "MLH1", "MSH2", 
                                   "MSH6", "PALB2", "PMS2"), 
                         cancers = c("BRA", "BC", "COL", "ENDO", 
                                     "GAS", "KID", "MELA", "OC", 
                                     "PANC", "PROS", "SMA"), 
                         markers = c("ER", "CK14", "CK5.6", 
                                     "PR", "HER2", "MSI"))
