library(abind)

K = 100
# Read in bootstrap diagnostics
all_diagnostics_boot = NULL
for (k in 1:K){
  load(paste0("results/diagnostics/diagnostics_boot", k, ".rData"))
  all_diagnostics_boot = c(all_diagnostics_boot, diagnostics_boot)
}

# Pull out and reshape relevant diagnostics
reshaped_boot_diagnostics = lapply(1:dim(all_diagnostics_boot[[1]])[3], function(i){
  abind(lapply(all_diagnostics_boot, function(x){ x[,,i] }), along=3)
})

# Get summary statistics of bootstrap diagnostics
boot_diagnostics_summary = lapply(reshaped_boot_diagnostics, function(x){
  abind(apply(x, c(1,2), function(y){summary(y)[1:6]}), 
        SD=array(apply(x, c(1,2), sd), dim=c(1,dim(reshaped_boot_diagnostics[[1]])[1:2])), 
        CI95lo=array(apply(x, c(1,2), quantile, probs=0.025, na.rm=TRUE), 
                     dim=c(1,dim(reshaped_boot_diagnostics[[1]])[1:2])), 
        CI95hi=array(apply(x, c(1,2), quantile, probs=0.975, na.rm=TRUE), 
                     dim=c(1,dim(reshaped_boot_diagnostics[[1]])[1:2])), 
        along=1)
})
names(boot_diagnostics_summary) = dimnames(all_diagnostics_boot[[1]])[[3]]

# Save and return summary statistics
save(boot_diagnostics_summary, 
     file="results/diagnostics/boot_diagnostics_summary.rData")
