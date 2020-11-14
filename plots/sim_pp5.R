library(abind)
library(knitr)

res_dir = "../sim/pp5/"
res_rm_dir = "../sim/pp5_riskmod/"
img_dir = "img/sim_pp5/"

# Read in and re-name diagnostics from models with risk modifiers
load(paste0(res_rm_dir, "results/diagnostics/diagnostics.rData"))
diagnostics_rm = diagnostics
load(paste0(res_rm_dir, "results/diagnostics/boot_diagnostics_summary.rData"))
boot_diagnostics_summary_rm = boot_diagnostics_summary

# Read in diagnostics
load(paste0(res_dir, "results/diagnostics/diagnostics.rData"))
load(paste0(res_dir, "results/diagnostics/boot_diagnostics_summary.rData"))

# Models
mods = c("PanelPRO"=1, "PanelPRO_sub"=2)

# Metrics
metrics = c("AUC"=1, "E/O"=2, "MSE"=5)
# Mutations
mutations = names(diagnostics)
names(mutations) = gsub("_", " ", mutations)


# Pull out and reshape relevant diagnostics
reshaped_diagnostics = lapply(metrics, function(x){
  rbind(sapply(mutations, function(y){ diagnostics[[y]][x,mods] }), 
        sapply(mutations, function(y){ diagnostics_rm[[y]][x,mods] }))[c(1,3,2,4),]
})

# Put all bootstrap confidence intervals together
# Lower CI
reshaped_diagnostics_CI95lo = 
  lapply(metrics, 
         function(x){
           rbind(sapply(mutations, function(y){
             boot_diagnostics_summary[[y]]["CI95lo",x,mods]
           }), 
           sapply(mutations, function(y){
             boot_diagnostics_summary_rm[[y]]["CI95lo",x,mods]
           }))[c(1,3,2,4),]
         })
# Upper CI
reshaped_diagnostics_CI95hi = 
  lapply(metrics, 
         function(x){
           rbind(sapply(mutations, function(y){
             boot_diagnostics_summary[[y]]["CI95hi",x,mods]
           }), 
           sapply(mutations, function(y){
             boot_diagnostics_summary_rm[[y]]["CI95hi",x,mods]
           }))[c(1,3,2,4),]
         })


# Colors for plotting
my_cols = c("#000000", "#0072B2")
my_cols_mat = matrix(rep(my_cols, each=length(mutations)*2), ncol=length(my_cols)*2)
# Point shapes for plotting
my_pch = c(16, 1)
my_pch_mat = matrix(rep(my_pch, times=length(mutations)*2), ncol=length(my_pch)*2, byrow=TRUE)
# Positions for mutations
y_mat = matrix(rep(1:length(mutations))+
                 rep(c(-0.3, -0.1, 0.1, 0.3), 
                     each=length(mutations)), ncol=length(my_cols)*2)

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
png(paste0(img_dir, "diagnostics.png"), width=600, height=260)
layout(mat = matrix(c(1, 2, 3, 
                      4, 4, 4), 
                    nrow = 2, ncol = 3, byrow = TRUE),
       heights = c(3.5, 1), widths = c(4, 3, 3))

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
  
  # y-axis mutation labels
  if (x=="AUC") {
    axis(2, 1:length(mutations), FALSE)
    text(y=1:length(mutations), x=par("usr")[1]-(xmax[x]-xmin[x])/20, 
         srt=50, adj=1, labels=names(mutations), xpd=TRUE, cex=1.5)
  }
  
  # Reference line for AUC and calibration
  if(x=="AUC" || x=="E/O" || x=="O/E") {
    segments(y0=0.65, y1=length(mutations)+0.35,
             x0=1, col="grey", lty=2)
  } else if (x=="E-O") {
    segments(y0=0.65, y1=length(mutations)+0.34,
             x0=0, col="grey", lty=2)
  }
}))

# Legends
plot.new()
par(bg=NA, xpd=TRUE)
legend("bottomleft", inset=c(0.2,-3), 
       legend=c("Full Model", "Submodel"), 
       fill=my_cols, bty="n", ncol=2, cex=1.4)
legend("bottomright", inset=c(0,-3), 
       legend=c("No Risk Modifiers", "Risk Modifiers"), 
       pch=my_pch, bty="n", ncol=2, lty=1, cex=1.4)
dev.off()


# Tables of values
lapply(names(metrics), function(x) {
  tab = data.frame(as.vector(sapply(names(mutations), function(y){ c(y, "", "", "")})),
                   c("Full model", "Full model", "Submodel", "Submodel"), 
                   c("No", "Yes"), 
                   as.vector(reshaped_diagnostics[[x]]), 
                   as.vector(reshaped_diagnostics_CI95lo[[x]]), 
                   as.vector(reshaped_diagnostics_CI95hi[[x]]))
  tab = tab[apply(tab, 1, function(x){!any(is.na(x))}),]
  names(tab) = c("Mutation", "Model", "Risk Modifiers", "Estimate", 
                 "Bootstrap 2.5%", "Bootstrap 97.5%")
  rownames(tab) = NULL
  return(kable(tab, format="latex", digits=5))
})
