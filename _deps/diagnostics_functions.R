# Discrimination (AUC)
library(pROC)

# Calibration (does it sum up to the number of carriers?)
# Ratio of the sum of predicted to the sum of observed
# Drops NAs
# x = vector of observed values
# pred = vector of predicted values
calibration = function(x, pred){
  na_idx = is.na(x) | is.na(pred)
  return((sum(pred[!na_idx])/sum(x[!na_idx])))
}

calibration_diff = function(x, pred){
  na_idx = is.na(x) | is.na(pred)
  return(sum(pred[!na_idx]) - sum(x[!na_idx]))
}

# Mean squared error
# Drops NAs
# x = vector of observed values
# pred = vector of predicted values
mse = function(x, pred) {
  return(mean((x-pred)^2, na.rm=TRUE))
}

# Obtain all three diagnostics
# x = vector of observed values
# pred = vector of predicted values
# title = character string to include in the plot title besides the AUC
# ... additional arguments to pass to roc.plot
get_diagnostics = function(x, pred) {
  return(c(AUC=tryCatch(auc(roc(x, pred)), error=function(err) NA),
           "calib. E/O"=calibration(x, pred),
           "calib. O/E"=calibration(pred, x),
           "calib. E-O"=calibration_diff(x, pred),
           MSE = mse(x, pred)))
}


get_all_diagnostics = function(gene, BM_model, 
                               mut_df, prob_df, noncarrier_idx, 
                               return_carriers=TRUE) {
  # Case 1: any gene
  if (length(gene)==1 && gene=="Any") {
    # Carriers of any mutation
    mutation_status = ifelse(rowSums(mut_df, na.rm = TRUE)==0, 0, 1)
    
    # Diagnostics with PanelPRO carrier probability of any mutation
    carrier_PP = colSums(prob_df["PanelPRO",-1,], na.rm = TRUE) 
    carrier_PP = carrier_PP / 
      (carrier_PP + prob_df["PanelPRO",1,])
    
    # Case 2: any BayesMendel gene
  } else if (length(gene)==1 && gene=="Any_BM") {
    # Vector of BayesMendel genes needs to be in the same order as 
    # the vector passed into PanelPRO
    gene = sort(intersect(colnames(prob_df), 
                          c("BRCA1", "BRCA2", "CDKN2A", 
                            "MLH1", "MSH2", "MSH6", "PANC")))
    if (is.null(BM_model)) {
      # BayesMendel models run
      BM_model = c("brcapro", "mmrpro", "melapro")
    }
    
    # Carriers of any gene in the vector
    mutation_status = ifelse(rowSums(mut_df[,gene], na.rm = TRUE)==0, 0, 1)
    
    # All genes and gene combinations
    all_genes = unique(sapply(gene, grep, colnames(prob_df)))
    
    # Diagnostics with PanelPRO carrier probability of 
    # any gene in the vector
    carrier_PP = colSums(prob_df["PanelPRO",all_genes,], na.rm = TRUE)
    carrier_PP = carrier_PP / 
      (carrier_PP + prob_df["PanelPRO",1,])
    
    # Diagnostics with PanelPRO version of BayesMendel model
    carrier_PP_sub = apply(prob_df[paste0("PP_", BM_model),all_genes,], 
                           3, sum, na.rm=TRUE)
    # Diagnostics with BayesMendel model
    carrier_BM_sub = apply(prob_df[paste0("BM_", BM_model),all_genes,], 
                           3, sum, na.rm=TRUE)
    
    # Case 3: pass in a single gene
  } else if (length(gene)==1) {
    # Carriers of the gene
    mutation_status = mut_df[,gene]
    
    # Diagnostics with PanelPRO carrier probability 
    # of single mutation for gene
    carrier_PP = prob_df["PanelPRO",gene,] / 
      (prob_df["PanelPRO",gene,] + prob_df["PanelPRO",1,])
    
    if (!is.null(BM_model)) {
      # Diagnostics with PanelPRO version of BayesMendel model
      carrier_PP_sub = prob_df[paste0("PP_", BM_model),gene,] / 
        (prob_df[paste0("PP_", BM_model),gene,] + 
           prob_df[paste0("PP_", BM_model),1,])
      # Diagnostics with BayesMendel model
      carrier_BM_sub = prob_df[paste0("BM_", BM_model),gene,] / 
        (prob_df[paste0("BM_", BM_model),gene,] + 
           prob_df[paste0("BM_", BM_model),1,])
    }
    
    # Case 4: pass in a vector of genes
  } else if (length(gene)>1) {
    # Carriers of any gene in the vector
    mutation_status = ifelse(rowSums(mut_df[,gene], na.rm = TRUE)==0, 0, 1)
    
    # All genes and gene combinations
    all_genes = unique(as.numeric(sapply(gene, grep, colnames(prob_df))))
    
    # Diagnostics with PanelPRO carrier probability of 
    # any gene in the vector
    carrier_PP = colSums(prob_df["PanelPRO",all_genes,], na.rm = TRUE)
    carrier_PP = carrier_PP / 
      (carrier_PP + prob_df["PanelPRO",1,])
    
    if (!is.null(BM_model)) {
      # Diagnostics with PanelPRO version of BayesMendel model
      carrier_PP_sub = colSums(prob_df[paste0("PP_", BM_model),-1,], na.rm=TRUE)
      # Diagnostics with BayesMendel model
      carrier_BM_sub = colSums(prob_df[paste0("BM_", BM_model),-1,], na.rm=TRUE)
    }
    
  } else {
    stop("No case works.")
  }
  
  # Indices of subjects to keep 
  keep_idx = noncarrier_idx | mutation_status==1
  
  # Diagnostics with PanelPRO model
  out_PP = get_diagnostics(mutation_status[keep_idx], carrier_PP[keep_idx])
  
  if (!is.null(BM_model)) {
    # Diagnostics with PanelPRO version of BayesMendel model
    out_PP_sub = get_diagnostics(mutation_status[keep_idx], 
                                 carrier_PP_sub[keep_idx])
    # Diagnostics with BayesMendel model
    out_BM_sub = get_diagnostics(mutation_status[keep_idx], 
                                 carrier_BM_sub[keep_idx])
    
    # Diagnostics ouput
    out = cbind(PanelPRO=out_PP, 
                PanelPRO_sub=out_PP_sub, 
                BayesMendel_sub=out_BM_sub)
    
    # Output with carriers
    if (return_carriers==TRUE) {
      out = list(diagnostics=out, 
                 carrier_PP=carrier_PP, 
                 carrier_PP_sub=carrier_PP_sub, 
                 carrier_BM_sub=carrier_BM_sub, 
                 mutation_status=mutation_status, 
                 keep_idx=keep_idx)
    } 
  } else {
    # Diagnostics ouput
    out = cbind(PanelPRO=out_PP, PanelPRO_sub=NA, BayesMendel_sub=NA)
    
    # Output with carriers
    if (return_carriers==TRUE) {
      out = list(diagnostics=out, 
                 carrier_PP=carrier_PP, 
                 mutation_status=mutation_status, 
                 keep_idx=keep_idx)
    } 
  }
  
  return(out)
}
