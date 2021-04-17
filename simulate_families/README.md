Functions for simulating detailed pedigrees, including family history of cancer, genotypes, and tumor marker testing results. The primary function, `sim.runSimFam`, can be found in `sim.runSimFam.R`; example usage is below. For further details, see the documentation for individual functions. 

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

# Simulate family using `PedUtils` code
fam = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                    PanelPRODatabase, genes, cancers, 
                    includeGeno = FALSE, includeBiomarkers = TRUE)
                    
# PanelPRO can be run on the simulated family
out = PanelPRO:::PanelPRO11(fam)
```
