library(abind)

# Read in diagnostics
load("results/diagnostics/diagnostics.rData")

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
        CI95lo=array(apply(x, c(1,2), function(y){
          quantile(y[is.finite(y)], probs=0.025, na.rm=TRUE)
        }), 
        dim=c(1,dim(reshaped_boot_diagnostics[[1]])[1:2])), 
        CI95hi=array(apply(x, c(1,2), function(y){
          quantile(y[is.finite(y)], probs=0.975, na.rm=TRUE)
        }), 
        dim=c(1,dim(reshaped_boot_diagnostics[[1]])[1:2])), 
        along=1)
})
names(boot_diagnostics_summary) = dimnames(all_diagnostics_boot[[1]])[[3]]

# Save and return summary statistics
save(boot_diagnostics_summary, 
     file="results/diagnostics/boot_diagnostics_summary.rData")



# Metrics
metrics = c("AUC"=1, "E/O"=2, "MSE"=5)
# Mutations
mutations = names(diagnostics)
names(mutations) = gsub("_", " ", mutations)


# Pull out and reshape relevant diagnostics
reshaped_diagnostics = lapply(metrics, function(x){
  sapply(mutations, function(y){ diagnostics[[y]][x,] })
})

# Put all bootstrap confidence intervals together
# Lower CI
reshaped_diagnostics_CI95lo = 
  lapply(metrics, 
         function(x){
           sapply(mutations, function(y){
             boot_diagnostics_summary[[y]]["CI95lo",x,]
           })
         })
# Upper CI
reshaped_diagnostics_CI95hi = 
  lapply(metrics, 
         function(x){
           sapply(mutations, function(y){
             boot_diagnostics_summary[[y]]["CI95hi",x,]
           })
         })


# Make plots comparing diagnostics
my_cols = c("black", "firebrick2", "dodgerblue2")
my_cols_mat = matrix(rep(my_cols, each=length(mutations)), ncol=3)
x_mat = matrix(rep(1:length(mutations), 3)+
                 rep(c(-0.15, 0, 0.15), 
                     each=length(mutations)), ncol=3)

# y-axis limits
finite.min = function(x){min(x[is.finite(x)], na.rm=TRUE)}
finite.max = function(x){max(x[is.finite(x)], na.rm=TRUE)}
ymin = sapply(c("AUC", "E/O", "MSE"), function(i){
  finite.min(as.numeric(c(reshaped_diagnostics_CI95lo[[i]], 
                          reshaped_diagnostics[[i]])))
})
ymax = sapply(c("AUC", "E/O", "MSE"), function(i){
  finite.max(as.numeric(c(reshaped_diagnostics_CI95hi[[i]], 
                          reshaped_diagnostics[[i]])))
})
ymax["AUC"] = 1


pdf("results/diagnostics/diagnostics.pdf", width=12, height=7)
par(mar=c(8.1, 4.1, 4.1, 2.1), xpd=TRUE)
invisible(sapply(names(metrics), function(x){
  plot(x_mat, t(reshaped_diagnostics[[x]]), 
       ylim=c(ymin[x], ymax[x]), xaxt="n",
       main=x, xlab="", ylab=x, 
       col=my_cols_mat, pch=16, 
       cex=1.5, cex.lab=1.3, cex.axis=1.3, cex.main=1.3)
  segments(x0=as.numeric(x_mat), x1=as.numeric(x_mat), 
           y0=as.numeric(t(reshaped_diagnostics_CI95lo[[x]])), 
           y1=as.numeric(t(reshaped_diagnostics_CI95hi[[x]])), 
           col=my_cols_mat, lwd=2)
  
  axis(1, 1:length(mutations), FALSE)
  text(x=1:length(mutations), y=par("usr")[3]-(ymax[x]-ymin[x])/20, 
       srt=30, adj=1, labels=names(mutations), xpd=TRUE, cex=1.3)
  
  # Reference line for AUC and calibration
  if(x=="AUC" || x=="E/O" || x=="O/E") {
    segments(x0=0.65, x1=length(mutations)+0.35,
             y0=1, col="grey", lty=2)
  } else if (x=="E-O") {
    segments(x0=0.65, x1=length(mutations)+0.34,
             y0=0, col="grey", lty=2)
  }
  
  # Legend
  legend("bottom", inset=c(0,-0.33), 
         legend=c("PanelPRO full", "PanelPRO sub", "BayesMendel sub"), 
         col=my_cols, pch=16, bty="n", ncol=3, lty=1)
}))
dev.off()


png("results/diagnostics/diagnostics%01d.png", width=650, height=330)
par(mar=c(7, 4.6, 2.8, 1), xpd=TRUE, bg=NA)
invisible(sapply(names(metrics), function(x){
  plot(x_mat, t(reshaped_diagnostics[[x]]), 
       ylim=c(ymin[x], ymax[x]), xaxt="n",
       main=x, xlab="", ylab=x, 
       col=my_cols_mat, pch=16, 
       cex=1.5, cex.lab=1.3, cex.axis=1.3, cex.main=1.3)
  segments(x0=as.numeric(x_mat), x1=as.numeric(x_mat), 
           y0=as.numeric(t(reshaped_diagnostics_CI95lo[[x]])), 
           y1=as.numeric(t(reshaped_diagnostics_CI95hi[[x]])), 
           col=my_cols_mat, lwd=2)
  
  axis(1, 1:length(mutations), FALSE)
  text(x=1:length(mutations), y=par("usr")[3]-(ymax[x]-ymin[x])/20, 
       srt=30, adj=1, labels=names(mutations), xpd=TRUE, cex=1.3)
  
  # Reference line for AUC and calibration
  if(x=="AUC" || x=="E/O" || x=="O/E") {
    segments(x0=0.65, x1=length(mutations)+0.35,
             y0=1, col="grey", lty=2)
  } else if (x=="E-O") {
    segments(x0=0.65, x1=length(mutations)+0.34,
             y0=0, col="grey", lty=2)
  }
  
  # Legend
  legend("bottom", inset=c(0,-0.53), 
         legend=c("PanelPRO full", "PanelPRO sub", "BayesMendel sub"), 
         col=my_cols, pch=16, bty="n", ncol=3, lty=1)
}))
dev.off()
