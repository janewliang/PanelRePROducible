res_dir = "../sim/pp11_extra/"
img_dir = "img/sim_pp11_extra/"

# Read in bootstrap diagnostics
load(paste0(res_dir, "results/diagnostics/boot_diagnostics_summary.rData"))

# Number of families
subset_n = c(1000000, 1200000, 1400000, 1600000, 1800000, 2000000)

# Genes
genes = names(boot_diagnostics_summary[[1]])
names(genes) = gsub("_", " ", genes)
# Genes that are included in the BayesMendel models
BM_genes = genes[c("BRCA1", "BRCA2", "BRCAPRO genes", 
                   "CDKN2A", "MLH1", "MSH2", "MSH6", 
                   "MMRpro genes", "Any BM")]
# Individual genes in the BayesMendel models (no groupings)
BM_genes_indiv = c("BRCA1", "BRCA2", "CDKN2A", 
                   "MLH1", "MSH2", "MSH6")

# Pull out and reshape bootstrap improvement frequencies
reshaped_if = sapply(improvement_freq, function(x) {
  sapply(BM_genes, 
         function(gene) {
           x[[gene]][1,"E/O"]
         })
})[,as.character(subset_n)]

# Colors for plotting
my_cols = c("#E69F00", "#009E73", "#882255", "#888888", 
            "#CC6677", "#56B4E9", "#F0E442", "#D55E00", 
            "#332288")

# Set up legend names
BM_genes_legend = names(BM_genes)
BM_genes_legend[c(3,8)] = c("BRCAPRO\ngenes", "MMRpro\ngenes")

# Plot percent of bootstraps that improve against number of simulated families
png(paste0(img_dir, "boot_if.png"), width=800, height=500)
par(mar=c(5.3, 5.3, 1, 8.3), bg=NA)
matplot(subset_n, t(reshaped_if), type = "b", 
        pch = 16, lty = 1, col = my_cols, 
        xlab = "Number of Families", 
        ylab = "Percent of Bootstraps that Improve", 
        cex=1.8, cex.lab=1.6, cex.axis=1.6, cex.main=1.8)

legend("topright", legend = BM_genes_legend, 
       fill=my_cols, title="Gene", bty="n", cex=1, 
       xpd=TRUE, inset=c(-0.17,0))
dev.off()


# Read in mutation counts and bootstraps
load(paste0(res_dir, "results/diagnostics/boot_mut_counts.rData"))
load(paste0(res_dir, "results/diagnostics/boot_diagnostics.rData"))

# Extract bootstraps for 2 million families
boot_2mill = reshaped_boot_diagnostics[[6]]
perc_diff_boot = sapply(BM_genes_indiv, function(gene){
                            (abs(boot_2mill[[gene]][2,2,]-1) - 
                               abs(boot_2mill[[gene]][2,1,]-1)) / 
                              (abs(boot_2mill[[gene]][2,2,]-1)+1)
                          })

# Plot bootstrap calibration against mutation carrier counts
png(paste0(img_dir, "calib_vs_carrier_counts_%01d.png"), 
    width = 350, height = 300)
par(mar = c(4, 4, 4, 1), bg=NA)
sapply(BM_genes_indiv, function(gene) {
           plot(boot_mut_counts[,gene], 
                boot_2mill[[gene]][2,1,], 
                cex = 0.8, pch = 16, col = "#0072B2", main = gene,
                xlab = "Number of Carriers", ylab = "E/O")
           points(boot_mut_counts[,gene], 
                  boot_2mill[[gene]][2,2,], 
                  cex = 0.8, pch = 16, col = "#E69F00")
           segments(x0 = boot_mut_counts[,gene], 
                    y0 = boot_2mill[[gene]][2,1,], 
                    y1 = boot_2mill[[gene]][2,2,], 
                    col = ifelse(perc_diff_boot[,gene] > 0, 
                                 "#0072B2", "#E69F00"))
           abline(h=1, col="grey")
         })
dev.off()

# Percent of bootstraps where PanelPRO has a larger calibration
sapply(BM_genes_indiv, function(gene) {
           mean(boot_2mill[[gene]][2,1,] >  
                  boot_2mill[[gene]][2,2,])
         })
