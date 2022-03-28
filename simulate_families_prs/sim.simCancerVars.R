#' Generate the \code{isAff} and \code{Age} Variables for Each Cancer
#'
#' Generates the \code{isAff} indicator and \code{Age} of diagnosis variables for 
#' each cancer, plus \code{isAffAny} (indicator for having any cancer) and 
#' \code{AgeAny} (age of diagnosis for first cancer). Also returns an overall 
#' \code{isDead} indicator variable and an updated \code{CurAge} variable based 
#' on the minimum of the death age and original current age. Incorporates 
#' censoring. Returns matrix with dimensions \code{number of members} by 
#' \code{2 x (number of cancers) + 4}. 
#' @param genoMat genotype matrix for the family (one row for each family member). 
#' @param Gender indicator vector (1 for males, 0 for females) of same length as 
#' number of rows in \code{genoMat}. 
#' @param CurAge vector of current ages for each family member; of same length as 
#' number of rows in \code{genoMat}. 
#' @param CP list of cancer penetrance matrices, one for each gender. 
#' @param PG vector of possible genotypes. 
#' @param cancers vector of cancers to include in the model. 
#' @param ageMax maximum age to consider. Defaults to 94. 
#' @param ageMin minimum age to consider. Defaults to 2. 
#' @param censoring if FALSE, then will assume everyone without any cancers lives
#' to the maximum age. Defaults to TRUE.
#' @param affTime if TRUE, the output will include the cancer times (which may 
#' or may not be observed due to censoring). Defaults to FALSE
#' @details Assumes naming conventions for cancers are consistent between arguments. 
#' The first column in the output matrix will be the current age \code{CurAge} and 
#' last column will be the indicator for death \code{isDead}. The 
#' \code{2 x (number of cancers) + 2} columns in the middle are the \code{isAff} 
#' and \code{Age} variables for each cancer, plus the two columns for any cancer. 
#' @family simulations
sim.simCancerVars = function(genoMat, Gender, CurAge, CP, PG, cancers, 
                             ageMax = 94, ageMin = 2, 
                             censoring = TRUE, affTime = FALSE,
                             latent = NULL, latentScore = NULL) {
  # Short cancer names
  cancers_short = PanelPRO:::CANCER_NAME_MAP$short[match(cancers, PanelPRO:::CANCER_NAME_MAP$long)]
  
  # Total number of people in family
  N = nrow(genoMat)
  
  if (length(Gender) != N) {
    stop("Gender should be of same length as number of people in genoMat")
  }
  if (length(CurAge) != N) {
    stop("CurAge should be of same length as number of people in genoMat")
  }
  
  # Break up the cancer density penetrance matrices
  cancerDens = list(CP$Dens["1",cancers,PG,,drop=FALSE], # First person is female (mother)
                    CP$Dens["2",cancers,PG,,drop=FALSE]) # Second person is male (her son)
  
  # Genotype combinations for all people in family 
  genoCombs = sapply(1:nrow(genoMat), function(i){
    paste(colnames(genoMat)[genoMat[i,] == 1], collapse=".")
  })
  genoCombs[genoCombs==""] = "noncarrier"
  # Figure out the index of each person's genotype combination compared to all possible
  genoCombsIdx = unlist(sapply(genoCombs, function(x){which(x == PG)}))
  
  if(censoring == TRUE){
    # Simulate age of death
    # Centered at 85, but values above ageMax are re-simulated to be centered around ageMax - 15. 
    deathAge = round(rnorm(N, 85, 10))
    while (any(deathAge > ageMax)) {
      deathAge[deathAge > ageMax] = round(rnorm(sum(deathAge > ageMax), ageMax - 15, 5))
    }
    
    # Set CurAge to the minimum of CurAge and deathAge 
    CurAge = pmin(CurAge, deathAge)
    
    # If deathAge is less than or equal to CurAge, the person is dead
    isDead = ifelse(deathAge <= CurAge, 1, 0)
  } else{
    isDead = rep(0, N)
  }
  
  # Set up matrix of isAffected variables, one for each cancer
  isAffectMat = matrix(nrow=N, ncol=length(cancers))
  colnames(isAffectMat) = paste0("isAff", cancers_short)
  # Set up matrix of Age variables, one for each cancer
  if(affTime == TRUE){
    ageMat = matrix(nrow=N, ncol=length(cancers)*2)
    colnames(ageMat) = c(paste0("Age", cancers_short), paste0("AffAge", cancers_short))
  } else{
    ageMat = matrix(nrow=N, ncol=length(cancers))
    colnames(ageMat) = paste0("Age", cancers_short)
  }
  
  for (j in 1:length(cancers)) { # Iterate through cancers
    # Times that each person got cancer j
    affectedTime = numeric(N) 
    
    for(i in 1:N){ # Iterate over each person
      # Simulate time person was affected by cancer
      pen = cancerDens[[Gender[i]+1]][,j,genoCombsIdx[i],,drop=FALSE] 
      pen[1:(length(pen)-1)] = pen[1:(length(pen)-1)] * latent$fact[latentScore[i]+1]
      pen[length(pen)] = max(1 - sum(pen[1:(length(pen)-1)]), 0)
      affectedTime[i] = max(which(rmultinom(1, 1, prob = pen)[,1] == 1), ageMin)
    }
    
    # If the person was affected at or before current age, the person is considered 
    # affected for that cancer
    isAffectMat[,j] = ifelse((affectedTime <= CurAge & 
                                affectedTime != ageMax + 1), 1, 0)
    # The age at which the person was affected should be the minimum of their current and 
    # affected ages (i.e. if person not affected, the age is just current age)

    ageMat[,j] = pmin(CurAge, affectedTime) 
    if(affTime == TRUE){
      ageMat[, j + length(cancers)] = affectedTime
    }
  }
  
  # Adding columns for having any cancer
  return(cbind(CurAge, 
               isAffectMat, isAffAny=apply(isAffectMat, 1, function(x){any(x==1)}), 
               ageMat, AgeAny=apply(ageMat, 1, min), isDead))
}
