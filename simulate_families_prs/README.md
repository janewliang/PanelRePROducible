Functions for simulating detailed pedigrees, including family history of cancer, genotypes, tumor marker testing results, and an additional PRS calculated for each individual based on an arbitrary number of independent SNPs. The primary function, `sim.runSimFam`, can be found in `sim.runSimFam.R`; example usage is below. For further details, see the documentation for individual functions. 

The function names are the same as those in the `simulate_families` directory, and their function is nearly identical aside from the additional PRS. 

```
# Cancers
cancers = c("Brain", "Breast", "Colorectal", "Endometrial", 
            "Gastric", "Kidney", "Melanoma", "Ovarian", 
            "Pancreas", "Prostate", "Small Intestine")
# Genes
genes = c("ATM", "BRCA1", "BRCA2", "CDKN2A", "CHEK2", "EPCAM", 
          "MLH1", "MSH2", "MSH6", "PALB2", "PMS2")
          
# Paternal aunts, paternal uncles
nSibsPatern = c(1, 2) 
# Maternal aunts, maternal uncles
nSibsMatern = c(2, 1) 
# Sisters and brothers
nSibs = c(2, 1) 
# We make the assumption that the number of sons and daughters for the 
# proband and all siblings, is the same. Nieces and nephews of the proband 
# are not sampled separately.
nGrandchild = c(1, 2) 

# Include a PRS calculated using 2 independent SNPs, each with allele frequency 0.1 
# The PRS score will be calculated as the sum of the SNP carrier statuses (0 for 
# noncarrier, 1 for carrier) and mapped to a multiplying factor that modifies 
# the cancer penetrances used to simulate each individual's cancer status, where 
# 0 corresponds to 0.8, 1 to 1, and 2 to 1.2. 
latent = data.frame(af = 0.1, 
                    fact = c(0.8, 1, 1.2))

# Simulate family using `PedUtils` code
fam = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                    PanelPRODatabase, genes, cancers, 
                    includeGeno = FALSE, includeBiomarkers = TRUE, latent = latent)
                    
# PanelPRO can be run on the simulated family
out = PanelPRO:::PanelPRO11(fam)
```
