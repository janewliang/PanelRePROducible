library(PanelPRO)
library(abind)

#' Underlying Function to Simulate a Family/Pedigree Matrix
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
#' @param alleleFreq allele frequencies for each gene of interest (named vector). 
#' @param CP list of cancer penetrance matrices, separated by each gender. 
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
#' @param BiomarkerTesting list of biomarker testing sensitivities and specificities.
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
#' @family simulations 
sim.simFam = function(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                      alleleFreq, CP, genes, cancers, 
                      maxMut = 2, ageMax = 94, ageMin = 2, 
                      includeGeno = FALSE, includeGrandparents = TRUE, 
                      censoring = TRUE, genderPro = NULL, 
                      genoMat = NULL, CurAge = NULL, affTime = FALSE, 
                      BiomarkerTesting = NULL, 
                      includeBiomarkers = FALSE, 
                      maxTries = 5) {
  # Check that specified cancers are in the CP object
  if (any(!(cancers %in% dimnames(CP$Dens)$cancers))) {
    stop("Cancers that are not in the CP object have been specified.")
  }
  
  # Check that specified genes are in the vector of allele frequencies
  if (any(!(genes %in% names(alleleFreq)))) {
    stop("Genes that are not in the vector of allele frequencies have been specified.")
  }
  # Subset allele frequencies to only include those genes
  alleleFreq = alleleFreq[genes]
  
  # Get possible genotypes and check that they are in the CP object
  PG = PanelPRO:::.getPossibleGenotype(genes, max_mut = maxMut)$list
  if (any(!(PG %in% dimnames(CP$Dens)$genotypes))) {
    stop("Genotypes that are not in the CP object have been specified.")
  }
  
  # Add the father to the number of male children in his branch
  if (length(nSibsPatern) != 2) {
    stop("nSibsPatern needs to be a numeric vector of length 2")
  }
  nChildPatern = nSibsPatern
  nChildPatern[2] = nChildPatern[2] + 1
  
  # Add the mother to the number of female children in her branch
  if (length(nSibsMatern) != 2) {
    stop("nSibsMatern needs to be a numeric vector of length 2")
  }
  nChildMatern = nSibsMatern
  nChildMatern[1] = nChildMatern[1] + 1
  
  # Add the proband to the number of female or male children in her branch/generation
  if (length(nSibs) != 2) {
    stop("nSibs needs to be a numeric vector of length 2")
  }
  nChild = nSibs
  if(is.null(genderPro)){
    genderPro <- ifelse(rbinom(1, 1, 0.5) == 1, "Female", "Male")
  }
  if(!is.null(genderPro) & !(genderPro %in% c("Female", "Male"))){
    genderPro <- ifelse(rbinom(1, 1, 0.5) == 1, "Female", "Male")
  }
  if(genderPro == "Female"){
    nChild[1] <- nChild[1] + 1
  } else{
    nChild[2] <- nChild[2] + 1
  }
  
  # If nGrandchild is passed in as a vector, use those numbers for each person in nChild
  if (is.null(nrow(nGrandchild))) {
    if (length(nGrandchild) != 2) {
      stop("If nGrandchild is passed in as a vector, it needs to have length 2")
    }
    nGrandchildInBranches = do.call(rbind, rep(list(nGrandchild), length=sum(nChild)))
    # Error if nGrandchild is a matrix with number of rows not equal to the number of 
    # people in nChild
  } else if (nrow(nGrandchild) != sum(nChild)) {
    stop("Number of rows in nGrandchild do not equal number of children in nChild.")
  } else {
    if (ncol(nGrandchild) != 2) {
      stop("If nGrandchild is passed in as a matrix, it needs to have 2 columns")
    }
    nGrandchildInBranches = nGrandchild
  }
  
  # Print warning message if includeGrandparents is FALSE but parents still have siblings
  if (includeGrandparents==FALSE && (nSibsPatern > 0 || nSibsMatern > 0)) {
    warning("You specified includeGrandparents==FALSE, but the mother and/or father 
            still has siblings. These will be included in the family pedigree as 
            unlinked members.")
  }
  
  
  # Generate genotype matrix for family
  if(is.null(genoMat)){
    genoMat = sim.buildGenoMat(alleleFreq, nChildPatern, nChildMatern, 
                               nChild, nGrandchildInBranches, maxTries)
  }
  # Flag to see if any family members have genotypes with more than 
  # maxMut mutations
  maxMutExceeded = any(rowSums(genoMat) > maxMut)
  # Intialize counter for re-runs
  counter = 1
  
  # Re-build genotype matrix until maxMut is not exceeded for anybody, 
  # or until the maximum number of runs has been hit. 
  while(maxMutExceeded == TRUE && counter < maxTries) {
    # Generate genotype matrix for family
    genoMat = sim.buildGenoMat(alleleFreq, nChildPatern, nChildMatern, 
                               nChild, nGrandchildInBranches, maxTries)
    
    # Flag to see if any family members have genotypes with more than 
    # maxMut mutations
    maxMutExceeded = any(rowSums(genoMat) > maxMut)
    # Increment counter
    counter = counter + 1
    # Print message
    print("Simulated genotype matrix contains genotypes with more than maxMut mutations, retrying.")
  }
  
  # Raise error if an admissible genotype matrix could not be simulated 
  # after maxTries attempts
  if (counter == maxTries) {
    stop(print(paste("Could not simulate an admissible genotype matrix for all family members after", 
                     maxTries, "attempts.")))
  } 
  
  # Total number of people in family
  N = nrow(genoMat)
  
  # Link Father ID, Mother ID, and genders for each branch of family
  # Paternal grandparents and their children 
  # (including father, who is the first male child)
  SexParentIDsPatern = sim.linkParents(nChildPatern[1], nChildPatern[2], 
                                          includeParents=includeGrandparents)
  lastID = nrow(SexParentIDsPatern) - 2*(!includeGrandparents)
  
  # Maternal grandparents and their children 
  # (including mother, who is the first female child)
  SexParentIDsMatern = sim.linkParents(nChildMatern[1], nChildMatern[2], 
                                          lastID=lastID, 
                                          includeParents=includeGrandparents)
  lastID = lastID + nrow(SexParentIDsMatern) - 2*(!includeGrandparents)
  
  # Proband and siblings (proband is the first female or male child listed)
  if(genderPro == "Female"){
    probandID <- lastID + 1
  } else{
    probandID <- lastID + nChild[1] + 1
  }
  SexParentIDsChild = sim.linkParents(nChild[1], nChild[2], 
                                         mothID=SexParentIDsMatern$ID[3], 
                                         fathID=SexParentIDsPatern$ID[3+nChildPatern[1]], 
                                         lastID=lastID)
  lastID = lastID + nrow(SexParentIDsChild)
  
  # Proband and siblings' spouses and children
  SexParentIDsGrandchild = data.frame()
  for (i in 1:nrow(nGrandchildInBranches)) {
    # Mother is blood relative in family, father is spouse
    if (SexParentIDsChild$Sex[i] == 0) {
      SexParentIDsGrandchild = rbind(SexParentIDsGrandchild, 
                                        sim.linkParents(nGrandchildInBranches[i,1], nGrandchildInBranches[i,2], 
                                                        mothID=SexParentIDsChild$ID[i], 
                                                        lastID=lastID+nrow(SexParentIDsGrandchild)))
      # Father is blood relative in family, mother is spouse
    } else {
      SexParentIDsGrandchild = rbind(SexParentIDsGrandchild, 
                                        sim.linkParents(nGrandchildInBranches[i,1], nGrandchildInBranches[i,2], 
                                                        fathID=SexParentIDsChild$ID[i], 
                                                        lastID=lastID+nrow(SexParentIDsGrandchild)))
    }
  }
  
  # Combine genders and parent IDs for all family members
  SexParentIDs = rbind(SexParentIDsPatern, SexParentIDsMatern,
                          SexParentIDsChild, SexParentIDsGrandchild)
  
  
  if(is.null(CurAge)){
    if(censoring == TRUE){
      # Generate current ages for each branch of family
      # Paternal grandparents and their children 
      # (including father, who is the first male child)
      CurAgePatern = sim.simCurAgeVar(sum(nChildPatern))
      
      # Maternal grandparents and their children 
      # (including mother, who is the first female child)
      CurAgeMatern = sim.simCurAgeVar(sum(nChildMatern))
      
      # Proband and siblings
      CurAgeChild = sim.simCurAgeVar(sum(nChild), CurAgePatern[3], CurAgeMatern[3])
      
      # Proband and siblings' spouses and children
      CurAgeGrandchild = unlist(lapply(1:nrow(nGrandchildInBranches), function(i){
        # Mother is blood relative in family, father is spouse
        if (SexParentIDsChild$Sex[i] == 0) {
          sim.simCurAgeVar(sum(nGrandchildInBranches[i,]), mothAge=CurAgeChild[i])
          # Father is blood relative in family, mother is spouse
        } else {
          sim.simCurAgeVar(sum(nGrandchildInBranches[i,]), fathAge=CurAgeChild[i])
        }
      }))
      
      # Combine current ages for all family members
      CurAge = c(CurAgePatern, CurAgeMatern, CurAgeChild, CurAgeGrandchild)
      
      # Ensure that all current ages fall between ageMin and ageMax
      CurAge[CurAge < ageMin] = ageMin
      CurAge[CurAge > ageMax] = ageMax
    } else{
      CurAge = rep(ageMax, N)
    }
  }
  
  # Get Affected and Age columns for each cancer, plus the isDead column
  Cancers = sim.simCancerVars(genoMat, SexParentIDs$Sex, CurAge, 
                              CP, PG, cancers, 
                              ageMax = ageMax, censoring = censoring,
                              affTime = affTime)
  
  
  
  # Designate the proband
  isProband = rep(0, length=N)
  isProband[probandID] = 1
  
  # Put the family pedigree together
  fam = data.frame(SexParentIDs, isProband, Cancers)
  
  # Drop the grandparents, if necessary
  if (includeGrandparents==FALSE) {
    fam = fam[fam$ID!=0,]
  }
  
  # Add other columns
  fam$Twins = 0
  fam$Ancestry = "nonAJ"
  fam$race = "All_Races"
  fam$interAge = fam$riskmod = list(character(0))
  
  
  # Simulate tumor markers
  if (includeBiomarkers == TRUE) {
    if (!(is.null(BiomarkerTesting))) {
      if ("Breast" %in% cancers) {
        if (all(c("BRCA1_hetero_anyPV", "BRCA2_hetero_anyPV") %in% genes)) {
          fam[, names(dim(BiomarkerTesting$Breast))[1:5]] = NA
          if (sum(fam$isAffBC == 1) > 0) {
            fam[fam$isAffBC == 1, names(dim(BiomarkerTesting$Breast))[1:5]] =
              t(apply(genoMat[fam$isAffBC == 1,c("BRCA1_hetero_anyPV", 
                                                 "BRCA2_hetero_anyPV"), 
                              drop = FALSE], 1, function(x) {
                probs = BiomarkerTesting$Breast[-1, -1, -1, -1, -1,
                                                as.character(x["BRCA1_hetero_anyPV"]),
                                                as.character(x["BRCA2_hetero_anyPV"])]
                # Need to rescale the probabilities according to how many marker result combos are included in the group
                probs = probs / apply(probs, 1:length(dim(probs)), function(x){sum(probs==x)})
                arrayInd(sample(1:length(probs), 1, replace = TRUE, prob = probs),
                         .dim=rep(2,5)) - 1
              }))
          }
        } else {
          warning("Need to include both BRCA1_hetero_anyPV and BRCA2_hetero_anyPV to simulate breast cancer markers.")
        }
      }  
      if ("Colorectal" %in% cancers) {
        if (all(c("MLH1_hetero_anyPV", "MSH2_hetero_anyPV", 
                  "MSH6_hetero_anyPV", "PMS2_hetero_anyPV") %in% genes)) {
          fam$MSI = NA
          if (sum(fam$isAffCOL == 1) > 0) {
            fam$MSI[fam$isAffCOL == 1] = apply(genoMat[fam$isAffCOL == 1,
                                                       c("MLH1_hetero_anyPV", 
                                                         "MSH2_hetero_anyPV", 
                                                         "MSH6_hetero_anyPV", 
                                                         "PMS2_hetero_anyPV"), 
                                                       drop = FALSE], 1, function(x) {
              probs = BiomarkerTesting$Colorectal[-1,
                                                  as.character(x["MLH1_hetero_anyPV"]),
                                                  as.character(x["MSH2_hetero_anyPV"]),
                                                  as.character(x["MSH6_hetero_anyPV"]),
                                                  as.character(x["PMS2_hetero_anyPV"])]
              sample(c(0,1), 1, replace = TRUE, prob = probs)
            })
          }
        }else {
          warning("Need to include all four of MLH1_hetero_anyPV, MSH2_hetero_anyPV, MSH6_hetero_anyPV and PMS2_hetero_anyPV to simulate colorectal cancer markers.")
        }
      } 
      
      if (length(intersect(c("Breast", "Colorectal"), cancers)) == 0) {
        warning("No cancers with biomarkers were specified, so skipping tumor marker simulations.")
      }
    } else {
      warning("No biomarker testing information supplied, so skipping tumor marker simulations.")
    }
  } 
  
  
  # Include the genotype matrix for all members, if necessary
  if (includeGeno==TRUE) {
    colnames(genoMat) = PanelPRO:::formatGeneNames(colnames(genoMat), 
                                                   format = "only_gene")
    fam = data.frame(fam, genoMat)
  } 
  
  return(fam)
}
