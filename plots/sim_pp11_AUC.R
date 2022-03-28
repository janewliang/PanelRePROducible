library(abind)
library(knitr)
library(RColorBrewer)

cur_dir = getwd()
res_dir = "../sim/pp11/"
img_dir = "img/sim_pp11/"

# Read in diagnostics
load(paste0(res_dir, "results/diagnostics/diagnostics.rData"))

# Load functions for running models
setwd(res_dir)
source("../../_deps/run_model_functions.R")
setwd(cur_dir)

###############################################################################

## Cancers and genes to use
# Short cancer names (including CBC)
cancers = c("BRA", "BC", "COL", "ENDO", "GAS", "KID", 
            "MELA", "OC", "PANC", "PROS", "SI", "CBC")
# Look up long cancer names (don't include CBC)
cancers_long = PanelPRO:::CANCER_NAME_MAP$long[sapply(cancers[-length(cancers)], function(x){
  which(x==PanelPRO:::CANCER_NAME_MAP$short)
})]
# Genes
genes = c("ATM", "BRCA1", "BRCA2", "CDKN2A", "CHEK2", "EPCAM", 
          "MLH1", "MSH2", "MSH6", "PALB2", "PMS2")


## Create a dummy family to generate penetrance densities and survivals
# Empty data frame of cancer affectation statuses
dummy.cancers = setNames(as.data.frame(matrix(0, nrow=2, ncol=length(cancers))), 
                         paste0("isAff", cancers))
# Empty data frame of cancer ages
dummy.ages = setNames(as.data.frame(matrix(NA, nrow=2, ncol=length(cancers))), 
                      paste0("Age", cancers))
# Dummy family
dummy.fam = data.frame(ID=c(1,2), 
                       MotherID=c(0,1), 
                       FatherID=c(0,0), 
                       Sex=c(0,1), 
                       isProband=c(0,1), 
                       Twins=c(0,0), 
                       Ancestry=rep("nonAJ", 2), 
                       CurAge=c(60,30), 
                       isDead=c(0,0), 
                       race=rep("All_Races", 2), 
                       dummy.cancers, 
                       dummy.ages)
# Assign no risk modifiers
dummy.fam$interAge = dummy.fam$riskmod = list(character(0))


## Get cancer penetrance densities and survivals from the dummy family
# Build a dummy database
dummy.db = buildDatabase(genes=genes, 
                         cancers=cancers_long, 
                         ppd=BackCompatibleDatabase)
dummy.db$Contralateral = BackCompatibleDatabase$Contralateral

# Run `checkFam` on the dummy family
dummy.fam.checked = checkFam(dummy.fam, dummy.db)$ped_list[[1]]

# Define possible genotypes
PGs = PanelPRO:::.getPossibleGenotype(
  dummy.db$MS$ALL_GENE_VARIANTS, max_mut = 2, 
  homo_genes = dummy.db$MS$HOMOZYGOUS)

# Extract genotypes with no more than 1 mutation
direct_fill_PGs <- unname(PGs$list[rowSums(PGs$df) < 2])
# Extract genotypes with 2 or more mutations
# Character strings used for naming
multi_PGs <- unname(PGs$list[rowSums(PGs$df) >= 2]) 
# In list format
multi_muts <- strsplit(PGs$list, split = "\\.")[rowSums(PGs$df) >= 2] 

# Cancer penetrance densities and survivals
CP = PanelPRO:::calcCancerPenetrance(
  dummy.fam.checked, dummy.fam.checked$ID, 
  dummy.db, sub_dens = NULL, 
  PGs,direct_fill_PGs, multi_PGs, multi_muts, 
  net=TRUE, consider.modification=FALSE)

# Extract allele frequencies from database
all_gene_variants = PanelPRO:::DEFAULT_VARIANTS[genes]
alleles = unique(sub("_.*_", "_", all_gene_variants))
alleleFreq = dummy.db$af[,"nonAJ"][alleles]
names(alleleFreq) = genes

# Age at which to obtain relative risks
age = 70
RR_female = sapply(cancers_long, function(canc) {
  CP$Dens["1", canc, dummy.db$MS$ALL_GENE_VARIANTS, as.character(age)] / 
    CP$Dens["1", canc, "noncarrier", as.character(age)]
})
dimnames(RR_female)[[1]] = genes
RR_male = sapply(cancers_long, function(canc) {
  CP$Dens["2", canc, dummy.db$MS$ALL_GENE_VARIANTS, as.character(age)] / 
    CP$Dens["2", canc, "noncarrier", as.character(age)]
})
dimnames(RR_male)[[1]] = genes


# Pull out AUC for each gene
AUCs = sapply(genes, function(x){
  diagnostics[[x]]["AUC", "PanelPRO"]
})


###############################################################################

# Colors for plotting
my_cols = brewer.pal(12, "Paired")[c(seq(1,9,by=2), seq(2,12,by=2))]
my_pch = c(0:6, 15:18)


# Plot AUC for each gene against allele frequencies 
text_x = alleleFreq + c(0, -0.00008, 0.00001, 0, -0.00005, 
                        0.00005, -0.00005, 0.0001, 0, 0, 0)
text_y = AUCs + c(0.02, 0.02, 0.02, 0.02, 0.02, -0.02, 
                  0.02, 0.02, 0.02, 0.02, 0.02)
png(paste0(img_dir, "AUC_prev.png"), width=700, height=500)
par(mar=c(5.3, 5.3, 1, 1), bg=NA)

# AUC vs prevalence
plot(alleleFreq, AUCs, col=my_cols, pch=8, 
     xlab="Pathogenic Variant Allele Frequency", ylab="AUC", 
     cex=1.8, cex.lab=1.6, cex.axis=1.6, cex.main=1.8)
text(text_x, text_y, names(alleleFreq), cex=1.4)
dev.off()


# AUC vs RR for all cancers and sexes
png(paste0(img_dir, "AUC_RR.png"), width=800, height=500)
par(mar=c(5.3, 5.3, 1, 8.3), bg=NA)
plot(1, type="n", xaxt = "n", 
     xlim=c(log10(min(RR_female, RR_male, na.rm=TRUE)), 
            log10(max(RR_female, RR_male, na.rm=TRUE))), 
     ylim=c(min(AUCs), 1), 
     xlab="Carrier/Noncarrier Risk of Cancer", ylab="AUC", 
     cex=1.8, cex.lab=1.6, cex.axis=1.6, cex.main=1.8)
invisible(lapply(1:length(cancers_long), function(i){
  points(log10(RR_female[genes,cancers_long[i]]), AUCs, 
         col=my_cols, pch=my_pch[i], cex=1.8)
  points(log10(RR_male[genes,cancers_long[i]]), AUCs, 
         col=my_cols, pch=my_pch[i], cex=1.8)
}))
axis(1, at = log10(c(1, 2, 4, 10, 20, 50, 100)), 
     labels = c(1, 2, 4, 10, 20, 50, 100), 
     cex=1.8, cex.lab=1.6, cex.axis=1.6)
legend("topright", legend=genes, 
       fill=my_cols, title="Gene", bty="n", cex=1, 
       xpd=TRUE, inset=c(-0.17,0))
legend("bottomright", legend=cancers_long, 
       pch=my_pch, title="Cancer", bty="n", cex=1, pt.cex=1.8, 
       xpd=TRUE, inset=c(-0.2,0))
dev.off()


# AUC vs RR for all cancers, female only
png(paste0(img_dir, "AUC_RR_female.png"), width=800, height=500)
par(mar=c(5.3, 5.3, 1, 8.3), bg=NA)
plot(1, type="n", xaxt = "n", 
     xlim=c(log10(min(RR_female, RR_male, na.rm=TRUE)), 
            log10(max(RR_female, RR_male, na.rm=TRUE))), 
     ylim=c(min(AUCs), 1), 
     xlab="Carrier/Noncarrier Risk of Cancer", ylab="AUC", 
     cex=1.8, cex.lab=1.6, cex.axis=1.6, cex.main=1.8)
invisible(lapply(1:length(cancers_long), function(i){
  points(log10(RR_female[genes,cancers_long[i]]), AUCs, 
         col=my_cols, pch=my_pch[i], cex=1.8)
}))
axis(1, at = log10(c(1, 2, 4, 10, 20, 50, 100)), 
     labels = c(1, 2, 4, 10, 20, 50, 100), 
     cex=1.8, cex.lab=1.6, cex.axis=1.6)
legend("topright", legend=genes, 
       fill=my_cols, title="Gene", bty="n", cex=1, 
       xpd=TRUE, inset=c(-0.17,0))
legend("bottomright", legend=cancers_long, 
       pch=my_pch, title="Cancer", bty="n", cex=1, pt.cex=1.8, 
       xpd=TRUE, inset=c(-0.2,0))
dev.off()


# AUC vs RR for all cancers, male only
png(paste0(img_dir, "AUC_RR_male.png"), width=800, height=500)
par(mar=c(5.3, 5.3, 1, 8.3), bg=NA)
plot(1, type="n", xaxt = "n", 
     xlim=c(log10(min(RR_female, RR_male, na.rm=TRUE)), 
            log10(max(RR_female, RR_male, na.rm=TRUE))), 
     ylim=c(min(AUCs), 1), 
     xlab="Carrier/Noncarrier Risk of Cancer", ylab="AUC", 
     cex=1.8, cex.lab=1.6, cex.axis=1.6, cex.main=1.8)
invisible(lapply(1:length(cancers_long), function(i){
  points(log10(RR_male[genes,cancers_long[i]]), AUCs, 
         col=my_cols, pch=my_pch[i], cex=1.8)
}))
axis(1, at = log10(c(1, 2, 4, 10, 20, 50, 100)), 
     labels = c(1, 2, 4, 10, 20, 50, 100), 
     cex=1.8, cex.lab=1.6, cex.axis=1.6)
legend("topright", legend=genes, 
       fill=my_cols, title="Gene", bty="n", cex=1, 
       xpd=TRUE, inset=c(-0.17,0))
legend("bottomright", legend=cancers_long, 
       pch=my_pch, title="Cancer", bty="n", cex=1, pt.cex=1.8, 
       xpd=TRUE, inset=c(-0.2,0))
dev.off()


# Tables of values
tab = data.frame(
  Estimate = c("AUC", "Allele Frequency", 
               as.vector(sapply(cancers_long, function(canc) { 
                 c(paste(canc, "RR", "(age 70)"), "")
               }))), 
  Sex = c("", "", rep(c("Female", "Male"), times=11)), 
  rbind(AUCs, alleleFreq, 
        do.call(rbind, lapply(cancers_long, function(canc) {
          rbind(RR_female[,canc], RR_male[,canc])
        }))))
tab = tab[apply(tab[-c(1,2)], 1, function(x){!any(is.nan(x))}),]
tab$Estimate[19] = "Prostate RR (age 70)"
rownames(tab) = NULL

kable(tab[,1:8], format="latex", digits=5)
kable(tab[,c(1:2, 9:13)], format="latex", digits=5)
