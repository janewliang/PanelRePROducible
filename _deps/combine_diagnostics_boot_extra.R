library(abind)

K = 1000
# Read in bootstrap diagnostics
subset_n = c(1000000, 1200000, 1400000, 1600000, 1800000, 2000000)
all_diagnostics_boot = lapply(subset_n, function(x) {
  vector(mode = "list", length = K)
})
names(all_diagnostics_boot) = subset_n
for (k in 1:K){
  load(paste0("results/diagnostics/diagnostics_boot", k, ".rData"))
  for (i in 1:length(all_diagnostics_boot)) {
    all_diagnostics_boot[[i]][[k]] = diagnostics_boot[[i]]
  }
}

# Pull out and reshape relevant diagnostics
reshaped_boot_diagnostics = lapply(all_diagnostics_boot, function(y) {
  reshaped_y = lapply(1:dim(y[[1]])[3], function(i){
    abind(lapply(y, function(x){ x[,,i] }), along=3)
  })
  names(reshaped_y) = dimnames(y[[1]])[[3]]
  return(reshaped_y)
})

# Get summary statistics of bootstrap diagnostics
boot_diagnostics_summary = lapply(reshaped_boot_diagnostics, function(z) {
  z_summary = lapply(z, function(x){
    abind(apply(x, c(1,2), function(y){summary(y)[1:6]}), 
          SD=array(apply(x, c(1,2), sd), dim=c(1,dim(z[[1]])[1:2])), 
          CI95lo=array(apply(x, c(1,2), quantile, probs=0.025, na.rm=TRUE), 
                       dim=c(1,dim(z[[1]])[1:2])), 
          CI95hi=array(apply(x, c(1,2), quantile, probs=0.975, na.rm=TRUE), 
                       dim=c(1,dim(z[[1]])[1:2])), 
          along=1)
  })
  names(z_summary) = dimnames(all_diagnostics_boot[[1]][[1]])[[3]]
  return(z_summary)
})

# Proportion of bootstraps where the full model improves over the sub-model
improvement_freq = lapply(reshaped_boot_diagnostics, function(y) {
  if_y = lapply(y, function(x) {
    rbind(PP = c("AUC" = mean(x[1,1,] >= x[1,2,]), 
                 "E/O" = mean(abs(1-x[2,1,]) <= abs(1-x[2,2,])), 
                 "MSE" = mean(x[5,1,] <= x[5,2,])), 
          BM = c("AUC" = mean(x[1,1,] >= x[1,3,]), 
                 "E/O" = mean(abs(1-x[2,1,]) <= abs(1-x[2,3,])), 
                 "MSE" = mean(x[5,1,] <= x[5,3,])))
  })
  names(if_y) = dimnames(all_diagnostics_boot[[1]][[1]])[[3]]
  return(if_y)
})

# Save bootstraps
save(reshaped_boot_diagnostics, 
     file="results/diagnostics/boot_diagnostics.rData")

# Save and return summary statistics
save(boot_diagnostics_summary, improvement_freq, 
     file="results/diagnostics/boot_diagnostics_summary.rData")
