library(PanelPRO)
library(abind)

#' Wrapper for Simulating a Family/Pedigree Matrix
#' 
#' Returns a completed family matrix with the following columns for each member: 
#' \itemize{
#'   \item \code{ID} = Member identifier
#'   \item \code{MotherID} = Mother's identifier number
#'   \item \code{FatherID} = Father's identifier number
#'   \item \code{Sex} = 1 for males, 0 for females
#'   \item \code{isProband} = Indicator for whether or not person is proband
#'   \item \code{CurAge} = Family member's current age, or age at death if dead
#'   \item \code{isAffCancer} variables = Indicators for whether or not person 
#'   had cancer type \code{Cancer}, e.g. \code{isAffBC}, \code{isAffO}, etc, plus
#'   an indicator \code{isAffAny} for having any cancer
#'   \item \code{AgeCancer} variables = Age of cancer diagnosis for cancer type 
#'   \code{Cancer}, e.g. \code{isAffBC}, \code{isAffO}, etc, plus the age of 
#'   first cancer diagnosis \code{AgeAny}. 0 if person did not have that cancer 
#'   type
#'   \item \code{isDead} = Indicator for whether or not person is dead
#'   \item \code{Twins} = Twin status, set to 0 for not twins (not simulated)
#'   \item \code{Ancestry} = Ancestry, set to "nonAJ" (not simulated)
#'   \item \code{race} = Race, set to "All_Races" (not simulated)
#'   \item \code{riskmod} = Risk modifier/intervention types, set to empty for 
#'   no risk modifiers (not simulated)
#'   \item \code{interAge} = Age of risk modifiers/interventions, set to 
#'   empty for no risk modifiers (not simulated)
#'   \item Optionally returns vectors for genotypes and tumor biomarkers, if 
#'   `includeGeno` and `includeBiomarkers` are set to `TRUE`. 
#' }
#' @param nSibsPatern vector of length \code{2}, indicating the number of sisters 
#' and brothers the father has (does not include father) 
#' @param nSibsMatern vector of length \code{2}, indicating the number of sisters 
#' and brothers the mother has (does not include mother) 
#' @param nSibs vector of length \code{2}, indicating the number of sisters and 
#' brothers the proband has (does not include proband) 
#' @param nGrandchild \code{sum(nSibs)+1} by \code{2} matrix, indicating the number 
#' of daughters and number of sons that the proband and each of her siblings have. 
#' Each row corresponds to the children of one of the members of \code{nChild}, with 
#' the proband as the first row. Alternatively, \code{nGrandchild} can be passed in 
#' as a vector of length \code{2}, indicating the number of daughters and number of 
#' sons that the proband and each of her siblings have, when they each have the same 
#' number.
#' @param database database to use, in the same format as `PanelPRODatabase` 
#' from the PanelPRO package
#' @param genes vector of genes to include in the model. 
#' @param cancers vector of cancers to include in the model. 
#' @param maxMut maximum number of mutations to consider. Defaults to 2. 
#' @param ageMax maximum age to consider. Defaults to 94.
#' @param ageMin minimum age to consider. Defaults to 2. 
#' @param includeGeno boolean flag, indicating whether or not to include genotype 
#' matrix in output. Defaults to FALSE, and was mostly used for troubleshooting. 
#' @param includeGrandparents boolean flag, indicating whether to drop the 
#' grandparents from the final pedigree. Defaults to TRUE. 
#' @param censoring if FALSE, then will assume everyone without any cancers lives
#' to the maximum age. Defaults to TRUE.
#' @param genderPro Can be set to "Female" or "Male" to designate the gender of the proband.
#' If left as NULL, the proband's gender will be randomly generated.
#' @param genoMat genotype matrix for the family. This is a matrix where the number of rows is
#' the number of family members and the number of columns is the number of mutations. If
#' left as NULL (the default), the genotypes will be generated.
#' @param CurAge vector of current ages for the family. If left as NULL (the default),
#' the current ages will be generated.
#' @param affTime if TRUE, the output will include the cancer times (which may 
#' or may not be observed due to censoring). Defaults to FALSE
#' @param includeBiomarkers boolean flag, indicating whether or not to simulate biomarkers 
#' (Breast and Colorectal currently supported). Defaults to FALSE. 
#' @param maxTries number of runs to attempt for simulating the allele 
#' array and genotype matrix for all family members. Pathogenic variants 
#' on both alleles for a given gene is currently not supported. When this 
#' occurs or an individual has a simulated genotype with more than 
#' maxMut simultaneous mutations, the function's solution is to re-simulate 
#' the alleles/genotypes up to maxTries times. If an admissible allele array 
#' or genotype matrix cannot be simulated after maxTries attempts, an 
#' error will be raised. These scenarios are more likely to occur for large 
#' values of alleleFreq. Defaults to 5 (each). 
#' @details Assumes naming conventions for cancers are consistent between arguments.
#' @family simulations export
#' @export
sim.runSimFam = function(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                         database, genes, cancers, 
                         maxMut = 2, ageMax = 94, ageMin = 2, 
                         includeGeno = FALSE, includeGrandparents = TRUE, 
                         censoring = TRUE, genderPro = NULL, 
                         genoMat = NULL, CurAge = NULL, affTime = FALSE, 
                         includeBiomarkers = FALSE, maxTries = 5, 
                         latent = NULL) {
  
  # Look up short cancer names and include CBC
  cancers_short = c(
    PanelPRO:::CANCER_NAME_MAP$short[sapply(cancers, function(x){
      which(x==PanelPRO:::CANCER_NAME_MAP$long)
    })], 
    "CBC"
  )
  
  
  ## Create a dummy family to generate penetrance densities and survivals
  # Empty data frame of cancer affectation statuses
  dummy.cancers = setNames(as.data.frame(matrix(0, nrow=2, ncol=length(cancers_short))), 
                           paste0("isAff", cancers_short))
  # Empty data frame of cancer ages
  dummy.ages = setNames(as.data.frame(matrix(NA, nrow=2, ncol=length(cancers_short))), 
                        paste0("Age", cancers_short))
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
                           cancers=cancers, 
                           ppd=database)
  dummy.db$Contralateral = database$Contralateral
  
  # Run `checkFam` on the dummy family
  dummy.fam.checked = checkFam(dummy.fam, dummy.db)$ped_list[[1]]
  
  # Define possible genotypes
  PGs = PanelPRO:::.getPossibleGenotype(
    dummy.db$MS$ALL_GENE_VARIANTS, max_mut = maxMut, 
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
  alleleFreq = dummy.db$af[alleles,"nonAJ"]
  names(alleleFreq) = all_gene_variants
  
  
  # Simulate family
  fam_PP = sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                      alleleFreq, CP, all_gene_variants, cancers, 
                      maxMut = maxMut, ageMax = ageMax, ageMin = ageMin, 
                      includeGeno = includeGeno, 
                      includeGrandparents = includeGrandparents, 
                      censoring = censoring, genderPro = genderPro, 
                      genoMat = genoMat, CurAge = CurAge, affTime = affTime, 
                      BiomarkerTesting = database$BiomarkerTesting, 
                      includeBiomarkers = includeBiomarkers, 
                      maxTries = maxTries, latent = latent)
}
