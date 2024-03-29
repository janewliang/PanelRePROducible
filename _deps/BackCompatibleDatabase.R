library(abind)
library(PanelPRO)
library(BayesMendel)

# New database
BackCompatibleDatabase = PanelPRODatabase


## Hypothetical PANC gene
# Add PANC allele frequency
BackCompatibleDatabase$AlleleFrequency = abind(BackCompatibleDatabase$AlleleFrequency, 
                                               PANC=rep(0.0034, dim(BackCompatibleDatabase$AlleleFrequency)[2]), 
                                               along=1)
names(dimnames(BackCompatibleDatabase$AlleleFrequency)) = c("Gene", "Ancestry")

# Set up empty PANC penetrance for each cancer
BackCompatibleDatabase$Penetrance = abind(BackCompatibleDatabase$Penetrance, 
                                          PANC=array(0, dim=dim(BackCompatibleDatabase$Penetrance[,1,,,,])), 
                                          along=2)
# Fill in PANC penetrances with SEER/non-carrier penetrances for each cancer
for (cancer in dimnames(BackCompatibleDatabase$Penetrance)$Cancer) {
  for (sex in c("Female", "Male")) {
    BackCompatibleDatabase$Penetrance[cancer,"PANC","All_Races",sex,,"Net"] = 
      BackCompatibleDatabase$Penetrance[cancer,"SEER","All_Races",sex,,"Net"]
  }
}
names(dimnames(BackCompatibleDatabase$Penetrance)) = c("Cancer", "Gene", "Race", 
                                                       "Sex", "Age", "PenetType")

# Fill in penetrances for developing pancreatic cancer conditional on PANC
BackCompatibleDatabase$Penetrance["Pancreas","PANC","All_Races","Female",,"Net"] = 
  penet.panc.net$fFX[,"P1"]
BackCompatibleDatabase$Penetrance["Pancreas","PANC","All_Races","Male",,"Net"] = 
  penet.panc.net$fMX[,"P1"]

# DOC
BackCompatibleDatabase$DOC = abind(BackCompatibleDatabase$DOC, 
                                   PANC=BackCompatibleDatabase$DOC[,"SEER",,,], 
                                   along=2)
names(dimnames(BackCompatibleDatabase$DOC)) = c("Cancer", "Gene", "Race", "Sex", "Age")

# Riskmod
BackCompatibleDatabase$Riskmod = abind(BackCompatibleDatabase$Riskmod, 
                                       PANC=array(1, dim=dim(BackCompatibleDatabase$Riskmod[,1,,,,])), 
                                       along=2)
names(dimnames(BackCompatibleDatabase$Riskmod)) = c("Cancer", "Gene", "Intervention", 
                                                    "Sex", "IntervAge", "DataType")

# Germline testing
BackCompatibleDatabase$GermlineTesting = abind(BackCompatibleDatabase$GermlineTesting, 
                                               PANC=rep(1, dim(BackCompatibleDatabase$GermlineTesting)[2]), 
                                               along=1)

## Set CBC penetrances to 0
dimsPP = dim(BackCompatibleDatabase$Contralateral)
BackCompatibleDatabase$Contralateral[1:dimsPP[1],1:dimsPP[2],1:dimsPP[3],1:dimsPP[4],1:dimsPP[5]] = 0




## Create a dummy family to generate penetrance densities and survivals
# Empty data frame of cancer affectation statuses
cancers11 = c("BRA", "BC", "COL", "ENDO", "GAS", "KID", "MELA", "OC", "PANC", "PROS", "SMA", "CBC")
cancers11_long = PanelPRO:::CANCER_NAME_MAP$long[sapply(cancers11[1:11], function(x){
  which(x==PanelPRO:::CANCER_NAME_MAP$short)
  })]
genes11 = c("ATM", "BRCA1", "BRCA2", "CDKN2A", "CHEK2", "EPCAM", "MLH1", "MSH2", "MSH6", "PALB2", "PMS2")
dummy.cancers = setNames(as.data.frame(matrix(0, nrow=4, ncol=12)), 
                         paste0("isAff", cancers11))
# Empty data frame of cancer ages
dummy.ages = setNames(as.data.frame(matrix(NA, nrow=4, ncol=12)), 
                      paste0("Age", cancers11))
# Dummy family
dummy.fam = data.frame(ID=c(1,2,3,4), 
                       MotherID=c(0,0,1,1), 
                       FatherID=c(0,0,2,2), 
                       Sex=c(0,1,0,1), 
                       isProband=c(0,0,0,1), 
                       Twins=c(0,0,0,0), 
                       Ancestry=rep("nonAJ", 4), 
                       CurAge=c(60,50,40,30), 
                       isDead=c(0,0,0,0), 
                       race=rep("All_Races", 4), 
                       dummy.cancers, 
                       dummy.ages)
dummy.fam$interAge = dummy.fam$riskmod = list(character(0))
# Give everybody breast cancer at various ages (for CBC)
dummy.fam$isAffBC = 1
dummy.fam$AgeBC = c(55, 45, 35, 25)
dummy.db = buildDatabase(genes=genes11, 
                         cancers=cancers11_long, 
                         ppd=BackCompatibleDatabase)
dummy.db$Contralateral = BackCompatibleDatabase$Contralateral
dummy.fam.checked = checkFam(dummy.fam, dummy.db)$ped_list[[1]]


# Get brcapro penetrances
# Database
brcapro.db = buildDatabase("BRCAPRO", ppd=BackCompatibleDatabase)

# Define possible genotypes
brcapro.PGs <- PanelPRO:::.getPossibleGenotype(
  brcapro.db$MS$ALL_GENE_VARIANTS,  max_mut = 2, 
  homo_genes = brcapro.db$MS$HOMOZYGOUS)

# Extract genotypes with no more than 1 mutation
brcapro.direct_fill_PGs <- unname(brcapro.PGs$list[rowSums(brcapro.PGs$df) < 2])
# Extract genotypes with 2 or more mutations
# Character strings used for naming
brcapro.multi_PGs <- unname(brcapro.PGs$list[rowSums(brcapro.PGs$df) >= 2]) 
# In list format
brcapro.multi_muts <- strsplit(brcapro.PGs$list, 
                               split = "\\.")[rowSums(brcapro.PGs$df) >= 2] 

# Cancer penetrance densities and survivals
brcapro.CP = PanelPRO:::calcCancerPenetrance(
  dummy.fam.checked, dummy.fam.checked$ID, 
  brcapro.db, sub_dens = NULL, 
  brcapro.PGs, brcapro.direct_fill_PGs, 
  brcapro.multi_PGs, brcapro.multi_muts, 
  net=TRUE, consider.modification=FALSE)


# mmrpro
# Database
mmrpro.db = buildDatabase("MMRPRO", ppd=BackCompatibleDatabase)

# Define possible genotypes
mmrpro.PGs <- PanelPRO:::.getPossibleGenotype(
  mmrpro.db$MS$ALL_GENE_VARIANTS,  max_mut = 2, 
  homo_genes = mmrpro.db$MS$HOMOZYGOUS)

# Extract genotypes with no more than 1 mutation
mmrpro.direct_fill_PGs <- unname(mmrpro.PGs$list[rowSums(mmrpro.PGs$df) < 2])
# Extract genotypes with 2 or more mutations
# Character strings used for naming
mmrpro.multi_PGs <- unname(mmrpro.PGs$list[rowSums(mmrpro.PGs$df) >= 2]) 
# In list format
mmrpro.multi_muts <- strsplit(mmrpro.PGs$list, 
                              split = "\\.")[rowSums(mmrpro.PGs$df) >= 2] 

# Cancer penetrance densities and survivals
mmrpro.CP = PanelPRO:::calcCancerPenetrance(
  dummy.fam.checked, dummy.fam.checked$ID, 
  mmrpro.db, sub_dens = NULL, 
  mmrpro.PGs, mmrpro.direct_fill_PGs, 
  mmrpro.multi_PGs, mmrpro.multi_muts, 
  net=TRUE, consider.modification=FALSE)


# pancpro
# Database
pancpro.db = buildDatabase(genes="PANC", cancers="Pancreas", 
                           ppd=BackCompatibleDatabase)

# Define possible genotypes
pancpro.PGs <- PanelPRO:::.getPossibleGenotype(
  pancpro.db$MS$ALL_GENE_VARIANTS,  max_mut = 2, 
  homo_genes = pancpro.db$MS$HOMOZYGOUS)

# Extract genotypes with no more than 1 mutation
pancpro.direct_fill_PGs <- unname(pancpro.PGs$list[rowSums(pancpro.PGs$df) < 2])
# Extract genotypes with 2 or more mutations
# Character strings used for naming
pancpro.multi_PGs <- unname(pancpro.PGs$list[rowSums(pancpro.PGs$df) >= 2]) 
# In list format
pancpro.multi_muts <- strsplit(pancpro.PGs$list, 
                               split = "\\.")[rowSums(pancpro.PGs$df) >= 2] 

# Cancer penetrance densities and survivals
pancpro.CP = PanelPRO:::calcCancerPenetrance(
  dummy.fam.checked, dummy.fam.checked$ID, 
  pancpro.db, sub_dens = NULL, 
  pancpro.PGs, pancpro.direct_fill_PGs, 
  pancpro.multi_PGs, pancpro.multi_muts, 
  net=TRUE, consider.modification=FALSE)


# melapro
# Database
melapro.db = buildDatabase(genes="CDKN2A", cancers="Melanoma", 
                           ppd=BackCompatibleDatabase)

# Define possible genotypes
melapro.PGs <- PanelPRO:::.getPossibleGenotype(
  melapro.db$MS$ALL_GENE_VARIANTS,  max_mut = 2, 
  homo_genes = melapro.db$MS$HOMOZYGOUS)

# Extract genotypes with no more than 1 mutation
melapro.direct_fill_PGs <- unname(melapro.PGs$list[rowSums(melapro.PGs$df) < 2])
# Extract genotypes with 2 or more mutations
# Character strings used for naming
melapro.multi_PGs <- unname(melapro.PGs$list[rowSums(melapro.PGs$df) >= 2]) 
# In list format
melapro.multi_muts <- strsplit(melapro.PGs$list, 
                               split = "\\.")[rowSums(melapro.PGs$df) >= 2] 

# Cancer penetrance densities and survivals
melapro.CP = PanelPRO:::calcCancerPenetrance(
  dummy.fam.checked, dummy.fam.checked$ID, 
  melapro.db, sub_dens = NULL, 
  melapro.PGs, melapro.direct_fill_PGs, 
  melapro.multi_PGs, melapro.multi_muts, 
  net=TRUE, consider.modification=FALSE)



# brcapro penetrances
penet.brca.net.PP = penet.brca.net
# brcapro genotypes
brca.genotypes.BM = c("B00", "B10", "B01", "B11")
# Corresponding PanelPRO genotypes
brca.genotypes.PP = dimnames(brcapro.CP$Dens)$genotypes

# First person is female
# Breast
penet.brca.net.PP$fFX[,brca.genotypes.BM] = t(brcapro.CP$Dens["1", "Breast",brca.genotypes.PP,-95])
penet.brca.net.PP$fFX[,!(colnames(penet.brca.net.PP$fFX) %in% brca.genotypes.BM)] = 0
# Ovarian
penet.brca.net.PP$fFY[,brca.genotypes.BM] = t(brcapro.CP$Dens["1", "Ovarian",brca.genotypes.PP,-95])
penet.brca.net.PP$fFY[,!(colnames(penet.brca.net.PP$fFY) %in% brca.genotypes.BM)] = 0
# Second person is male
# Breast
penet.brca.net.PP$fMX[,brca.genotypes.BM] = t(brcapro.CP$Dens["2", "Breast",brca.genotypes.PP,-95])
penet.brca.net.PP$fMX[,!(colnames(penet.brca.net.PP$fMX) %in% brca.genotypes.BM)] = 0
# Ovarian
penet.brca.net.PP$fMY[,brca.genotypes.BM] = t(brcapro.CP$Dens["2", "Ovarian",brca.genotypes.PP,-95])
penet.brca.net.PP$fMY[,!(colnames(penet.brca.net.PP$fMY) %in% brca.genotypes.BM)] = 0


# mmrpro penetrances
penet.mmr.net.PP = penet.mmr.net
# mmrpro genotypes
mmr.genotypes.BM = c("M000", "M100", "M010", "M001", "M110", "M101", "M011")
# Corresponding PanelPRO genotypes
mmr.genotypes.PP = dimnames(mmrpro.CP$Dens)$genotypes

# First person is female
# Colorectal
penet.mmr.net.PP$fFX[,mmr.genotypes.BM] = t(mmrpro.CP$Dens["1", "Colorectal",mmr.genotypes.PP,-95])
penet.mmr.net.PP$fFX[,!(colnames(penet.mmr.net.PP$fFX) %in% mmr.genotypes.BM)] = 0
# Endometrial
penet.mmr.net.PP$fFY[,mmr.genotypes.BM] = t(mmrpro.CP$Dens["1", "Endometrial",mmr.genotypes.PP,-95])
penet.mmr.net.PP$fFY[,!(colnames(penet.mmr.net.PP$fFY) %in% mmr.genotypes.BM)] = 0
# Second person is male
# Colorectal
penet.mmr.net.PP$fMX[,mmr.genotypes.BM] = t(mmrpro.CP$Dens["2", "Colorectal",mmr.genotypes.PP,-95])
penet.mmr.net.PP$fMX[,!(colnames(penet.mmr.net.PP$fMX) %in% mmr.genotypes.BM)] = 0
# Endometrial
penet.mmr.net.PP$fMY[,mmr.genotypes.BM] = t(mmrpro.CP$Dens["2", "Endometrial",mmr.genotypes.PP,-95])
penet.mmr.net.PP$fMY[,!(colnames(penet.mmr.net.PP$fMY) %in% mmr.genotypes.BM)] = 0


# pancpro penetrances
penet.panc.net.PP = penet.panc.net
# pancpro genotypes
panc.genotypes.BM = c("P0", "P1")
# Corresponding PanelPRO genotypes
panc.genotypes.PP = dimnames(pancpro.CP$Dens)$genotypes

# First person is female
# Pancreatic
penet.panc.net.PP$fFX[,panc.genotypes.BM] = t(pancpro.CP$Dens["1", "Pancreas",panc.genotypes.PP,-95])
penet.panc.net.PP$fFX[,!(colnames(penet.panc.net.PP$fFX) %in% panc.genotypes.BM)] = 0
# Second person is male
# Pancreatic
penet.panc.net.PP$fMX[,panc.genotypes.BM] = t(pancpro.CP$Dens["2", "Pancreas",panc.genotypes.PP,-95])
penet.panc.net.PP$fMX[,!(colnames(penet.panc.net.PP$fMX) %in% panc.genotypes.BM)] = 0


# melapro penetrances
penet.mela.hbi.net.PP = penet.mela.hbi.net
# melapro genotypes
mela.genotypes.BM = c("P160", "P161")
# Corresponding PanelPRO genotypes
mela.genotypes.PP = dimnames(melapro.CP$Dens)$genotypes

# First person is female
# Melanoma
penet.mela.hbi.net.PP$fFX[,mela.genotypes.BM] = t(melapro.CP$Dens["1", "Melanoma",mela.genotypes.PP,-95])
penet.mela.hbi.net.PP$fFX[,!(colnames(penet.mela.hbi.net.PP$fFX) %in% mela.genotypes.BM)] = 0
# Second person is male
# Melanoma
penet.mela.hbi.net.PP$fMX[,mela.genotypes.BM] = t(melapro.CP$Dens["2", "Melanoma",mela.genotypes.PP,-95])
penet.mela.hbi.net.PP$fMX[,!(colnames(penet.mela.hbi.net.PP$fMX) %in% mela.genotypes.BM)] = 0


## Set brcapro CBC to be the same as PanelPRODatabase
# # Contralateral breast cancer penetrances
# CBRCApenet.2016.PP = CBRCApenet.2016
# # brcapro genotypes
# cbc.genotypes.BM = c("G000", "G100", "G010", "G110")
# # Corresponding PanelPRO genotypes
# cbc.genotypes.PP = c("noncarrier", "BRCA1", "BRCA2", "BRCA1.BRCA2")
# 
# # First person is a female over 40
# CBRCApenet.2016.PP$fFX.Over40[1:95, cbc.genotypes.BM] = t(brcapro.CP$Dens["1", "Contralateral", cbc.genotypes.PP,])
# CBRCApenet.2016.PP$fFX.Over40[96:111, cbc.genotypes.BM] = 0
# CBRCApenet.2016.PP$fFX.Over40[,!(colnames(CBRCApenet.2016.PP$fFX.Over40) %in% cbc.genotypes.BM)] = 0
# # First person is a male over 40
# CBRCApenet.2016.PP$fMX.Over40[1:95, cbc.genotypes.BM] = t(brcapro.CP$Dens["2", "Contralateral", cbc.genotypes.PP,])
# CBRCApenet.2016.PP$fMX.Over40[96:111, cbc.genotypes.BM] = 0
# CBRCApenet.2016.PP$fMX.Over40[,!(colnames(CBRCApenet.2016.PP$fMX.Over40) %in% cbc.genotypes.BM)] = 0
# # First person is a female under 40
# CBRCApenet.2016.PP$fFX.Under40[1:95, cbc.genotypes.BM] = t(brcapro.CP$Dens["3", "Contralateral", cbc.genotypes.PP,])
# CBRCApenet.2016.PP$fFX.Under40[96:111, cbc.genotypes.BM] = 0
# CBRCApenet.2016.PP$fFX.Under40[,!(colnames(CBRCApenet.2016.PP$fFX.Under40) %in% cbc.genotypes.BM)] = 0
# # First person is a male under 40
# CBRCApenet.2016.PP$fMX.Under40[1:95, cbc.genotypes.BM] = t(brcapro.CP$Dens["4", "Contralateral", cbc.genotypes.PP,])
# CBRCApenet.2016.PP$fMX.Under40[96:111, cbc.genotypes.BM] = 0
# CBRCApenet.2016.PP$fMX.Under40[,!(colnames(CBRCApenet.2016.PP$fMX.Under40) %in% cbc.genotypes.BM)] = 0



## Set CBC to 0 for both brcapro and PanelPRO
# Contralateral breast cancer penetrances
cbc.net.PP = cbc.net
dimsBM = dim(cbc.net.PP$fFX.Over40)

# First person is a female over 40
cbc.net.PP$fFX.Over40[1:dimsBM[1],1:dimsBM[2]] = 0
# First person is a male over 40
cbc.net.PP$fMX.Over40[1:dimsBM[1],1:dimsBM[2]] = 0
# First person is a female under 40
cbc.net.PP$fFX.Under40[1:dimsBM[1],1:dimsBM[2]] = 0
# First person is a male under 40
cbc.net.PP$fMX.Under40[1:dimsBM[1],1:dimsBM[2]] = 0



# Save output
save(BackCompatibleDatabase, 
     penet.brca.net.PP, penet.mmr.net.PP, 
     penet.panc.net.PP, penet.mela.hbi.net.PP, 
     cbc.net.PP, 
     file="BackCompatibleDatabase.rData")
