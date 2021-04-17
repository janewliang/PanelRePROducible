library(abind)
library(knitr)

img_dir = "img/upper_pp11/"

# List of all directories to search in
all_dirs = c("../usc/pp11/", "../sim/pp11/", "../big_fam/pp11/")

# Read in diagnostics from each directory
diagnostics_list = lapply(all_dirs, function(dir_name){
  temp = try(load(paste0(dir_name, "results/diagnostics/subset_diagnostics.rData")), 
             silent=TRUE)
  if (class(temp) == "try-error") {
    load(paste0(dir_name, "results/diagnostics/diagnostics.rData"))
  }
  return(diagnostics)
})

# Read in summary of bootstrap diagnostics from each directory
boot_diagnostics_list = lapply(all_dirs, function(dir_name){
  temp = try(load(paste0(dir_name, "results/diagnostics/subset_boot_diagnostics_summary.rData")), 
             silent=TRUE)
  if (class(temp) == "try-error") {
    load(paste0(dir_name, "results/diagnostics/boot_diagnostics_summary.rData"))
  }
  return(boot_diagnostics_summary)
})

# Models
mods = c("PanelPRO"=1)

# Metrics
metrics = c("AUC"=1, "E/O"=2, "MSE"=5)
# Genes
genes = names(diagnostics_list[[1]])
names(genes) = gsub("_", " ", genes)


# Pull out and reshape relevant diagnostics
reshaped_diagnostics = lapply(metrics, function(x){
  rbind(sapply(genes, function(y){ diagnostics_list[[1]][[y]][x,mods] }), 
        sapply(genes, function(y){ diagnostics_list[[2]][[y]][x,mods] }), 
        sapply(genes, function(y){ diagnostics_list[[3]][[y]][x,mods] }))
})

# Put all bootstrap confidence intervals together
# Lower CI
reshaped_diagnostics_CI95lo = 
  lapply(metrics, 
         function(x){
           rbind(sapply(genes, function(y){
             boot_diagnostics_list[[1]][[y]]["CI95lo",x,mods]
           }), 
           sapply(genes, function(y){
             boot_diagnostics_list[[2]][[y]]["CI95lo",x,mods]
           }), 
           sapply(genes, function(y){
             boot_diagnostics_list[[3]][[y]]["CI95lo",x,mods]
           }))
         })
# Upper CI
reshaped_diagnostics_CI95hi = 
  lapply(metrics, 
         function(x){
           rbind(sapply(genes, function(y){
             boot_diagnostics_list[[1]][[y]]["CI95hi",x,mods]
           }), 
           sapply(genes, function(y){
             boot_diagnostics_list[[2]][[y]]["CI95hi",x,mods]
           }), 
           sapply(genes, function(y){
             boot_diagnostics_list[[3]][[y]]["CI95hi",x,mods]
           }))
         })


# Colors for plotting
my_cols = c("#000000", "#009E73", "#CC79A7")
my_cols_mat = matrix(rep(my_cols, each=length(genes)), ncol=length(my_cols))
# Point shapes for plotting
my_pch = c(16, 2, 4)
my_pch_mat = matrix(rep(my_pch, times=length(genes)), ncol=length(my_pch), byrow=TRUE)
# Positions for genes
y_mat = matrix(rep(1:length(genes))+
                 rep(c(-0.2, 0, 0.2), 
                     each=length(genes)), ncol=length(my_cols))

# x-axis limits
finite.min = function(x){min(x[is.finite(x)], na.rm=TRUE)}
finite.max = function(x){max(x[is.finite(x)], na.rm=TRUE)}
xmin = sapply(c("AUC", "E/O", "MSE"), function(i){
  finite.min(as.numeric(c(reshaped_diagnostics_CI95lo[[i]], 
                          reshaped_diagnostics[[i]])))
})
xmax = sapply(c("AUC", "E/O", "MSE"), function(i){
  finite.max(as.numeric(c(reshaped_diagnostics_CI95hi[[i]], 
                          reshaped_diagnostics[[i]])))
})
xmax["AUC"] = 1
# xlab labels/subfigure captions
xlab_name = c("AUC"="(a) Area under the curve (AUC)", 
              "E/O"="(b) Calibration (expected/observed)", 
              "MSE"="(c) Mean squared error (MSE)")


# Forest plots
png(paste0(img_dir, "diagnostics.png"), width=600, height=420)
layout(mat = matrix(c(1, 2, 3, 
                      4, 4, 4), 
                    nrow = 2, ncol = 3, byrow = TRUE),
       heights = c(7, 1), widths = c(4, 3, 3))

invisible(sapply(names(metrics), function(x){
  # Set up margins
  if (x=="AUC") {
    par(mar = c(4, 6.7, 1, 1), bg=NA)
  } else {
    par(mar = c(4, 0, 1, 1), bg=NA)
  }
  
  # Plot points and CIs
  plot(t(reshaped_diagnostics[[x]]), y_mat, 
       xlim=c(xmin[x], xmax[x]), yaxt="n",
       xlab=xlab_name[x], ylab="", 
       col=my_cols_mat, pch=my_pch_mat, 
       cex=1.5, cex.lab=1.4, cex.axis=1.5)
  segments(x0=as.numeric(t(reshaped_diagnostics_CI95lo[[x]])), 
           x1=as.numeric(t(reshaped_diagnostics_CI95hi[[x]])), 
           y0=as.numeric(y_mat), y1=as.numeric(y_mat), 
           col=my_cols_mat)
  
  # y-axis gene labels
  if (x=="AUC") {
    axis(2, 1:length(genes), FALSE)
    text(y=1:length(genes), x=par("usr")[1]-(xmax[x]-xmin[x])/20, 
         srt=60, adj=1, labels=names(genes), xpd=TRUE, cex=1.5)
  }
  
  # Reference line for AUC and calibration
  if(x=="AUC" || x=="E/O" || x=="O/E") {
    segments(y0=0.65, y1=length(genes)+0.35,
             x0=1, col="grey", lty=2)
  } else if (x=="E-O") {
    segments(y0=0.65, y1=length(genes)+0.34,
             x0=0, col="grey", lty=2)
  }
}))

# Legend
plot.new()
par(bg=NA, xpd=TRUE)
legend("bottom", inset=c(0,-6), 
       legend=c("HCP Cohort", "Simulations", "Upper Bound Simulations"), 
       col=my_cols, pch=my_pch, 
       bty="n", ncol=3, lty=1, cex=1.4)
dev.off()


# Tables of diagnostic metrics
lapply(names(metrics), function(x) {
  tab = data.frame(as.vector(sapply(names(genes), function(y){ c(y, "", "")})),
                   c("HCP Cohort", "Simulations", "Upper Bound Simulations"), 
                   as.vector(reshaped_diagnostics[[x]]), 
                   as.vector(reshaped_diagnostics_CI95lo[[x]]), 
                   as.vector(reshaped_diagnostics_CI95hi[[x]]))
  names(tab) = c("Gene", "Data", "Estimate", 
                 "Bootstrap 2.5%", "Bootstrap 97.5%")
  return(kable(tab, format="latex", digits=5))
})
