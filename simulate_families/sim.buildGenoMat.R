#' Build the Allele Array for All Family Members in a 4-Generation Family
#' 
#' Returns a 3D array of alleles, with dimensions \code{number of members} 
#' by \code{2} by \code{number of genes}
#' @param alleleFreq allele frequencies for each gene of interest
#' @param nChildPatern vector of length \code{2}, indicating the number of daughters 
#' (father's sisters) and number of sons (father + his brothers) in the paternal branch
#' @param nChildMatern vector of length \code{2}, indicating the number of daughters 
#' (mother + her sisters) and number of sons (mother's brothers) in the maternal branch
#' @param nChild vector of length \code{2}, indicating the number of daughters 
#' (proband + her sisters) and number of sons (proband's brothers) in the proband's 
#' generation
#' @param nGrandchildInBranches \code{2} by \code{sum(nChild)} matrix, indicating 
#' the number of daughters and number of sons that the proband and each of her siblings 
#' have. Each row corresponds to the children of one of the members of \code{nChild}, 
#' with the proband as the first row.
#' @family simulations
#' @importFrom abind abind
sim.buildAlleleArray = function(alleleFreq, nChildPatern, nChildMatern, 
                            nChild, nGrandchildInBranches){ 
  # Simulate allele matrices for paternal grandparents and their children 
  # (including father, who is the first male child)
  allelesPatern = sim.buildBranchOfAlleleMats(alleleFreq, sum(nChildPatern))
  
  # Simulate allele matrices for maternal grandparents and their children 
  # (including mother, who is the first female child)
  allelesMatern = sim.buildBranchOfAlleleMats(alleleFreq, sum(nChildMatern))
  
  # Simulate allele matrices for proband and siblings using parents' allele matrices
  # Proband is the first child and is female
  allelesChild = sim.buildBranchOfAlleleMats(alleleFreq, sum(nChild),
                                             matrix(allelesMatern[3,,], nrow=2), 
                                             matrix(allelesPatern[3+nChildPatern[1],,], nrow=2))
  
  # Simulate allele matrices for each of proband and siblings' spouses and children
  # Proband's children come first
  allelesGrandchildInBranches = abind(lapply(1:nrow(nGrandchildInBranches), function(i){
    sim.buildBranchOfAlleleMats(alleleFreq, sum(nGrandchildInBranches[i,]), matrix(allelesChild[i,,], nrow=2))
  }), along=1)
  
  # Put everything into the same 3d matrix
  alleleArray = abind(allelesPatern, allelesMatern, allelesChild, 
                      allelesGrandchildInBranches, along=1)
  
  return(alleleArray)
}


#' Build the Genotype Matrix for All Family Members in a 4-Generation Family
#' 
#' Returns a matrix of genotypes, with dimensions \code{number of members} 
#' by \code{number of genes}
#' @param alleleFreq allele frequencies for each gene of interest
#' @param nChildPatern vector of length \code{2}, indicating the number of daughters 
#' (father's sisters) and number of sons (father + his brothers) in the paternal branch
#' @param nChildMatern vector of length \code{2}, indicating the number of daughters 
#' (mother + her sisters) and number of sons (mother's brothers) in the maternal branch
#' @param nChild vector of length \code{2}, indicating the number of daughters 
#' (proband + her sisters) and number of sons (proband's brothers) in the proband's 
#' generation
#' @param nGrandchildInBranches \code{2} by \code{sum(nChild)} matrix, indicating 
#' the number of daughters and number of sons that the proband and each of her siblings 
#' have. Each row corresponds to the children of one of the members of \code{nChild}, 
#' with the proband as the first row.
#' @param maxTries number of runs to attempt for simulating the allele 
#' array and genotype matrix for all family members. Pathogenic variants 
#' on both alleles for a given gene is currently not supported. When this 
#' occurs or an individual has a simulated genotype with more than 
#' maxMut simultaneous mutations, the function's solution is to re-simulate 
#' the alleles/genotypes up to maxTries times. If an admissible allele array 
#' or genotype matrix cannot be simulated after maxTries attempts, an 
#' error will be raised. These scenarios are more likely to occur for large 
#' values of alleleFreq. Defaults to 5 (each). 
#' @family simulations
#' @importFrom abind abind
sim.buildGenoMat = function(alleleFreq, nChildPatern, nChildMatern, 
                            nChild, nGrandchildInBranches, maxTries){ 
  if (length(nChildPatern) != 2) {
    stop("nChildPatern needs to be a numeric vector of length 2")
  }
  if (length(nChildMatern) != 2) {
    stop("nChildMatern needs to be a numeric vector of length 2")
  }
  if (length(nChild) != 2) {
    stop("nChild needs to be a numeric vector of length 2")
  }
  if (ncol(nGrandchildInBranches) != 2) {
    stop("nGrandchildInBranches needs to have 2 columns")
  }
  if (nrow(nGrandchildInBranches) != sum(nChild)) {
    stop("Number of rows in nGrandchildInBranches does not equal number of children in nChild")
  }
  
  # Simulate 3D array of alleles
  alleleArray = sim.buildAlleleArray(alleleFreq, nChildPatern, 
                                     nChildMatern, nChild, 
                                     nGrandchildInBranches)
  # Sum array of alleles across the alleles
  alleleSumMat = apply(alleleArray, 1, function(x) {
    apply(x, 2, sum)
  })
  
  # Flag to see if any family members have pathogenic variants on both 
  # alleles for a given gene
  containsHomo = any(colSums(alleleSumMat) > 2)
  # Intialize counter for re-runs
  counter = 1
  
  # Re-build allele array until nobody has any homozygous mutations 
  # (currently not supported/considered non-viable), or until the 
  # maximum number of runs has been hit. 
  while(containsHomo == TRUE && counter < maxTries) {
    # Simulate 3D array of alleles
    alleleArray = sim.buildAlleleArray(alleleFreq, nChildPatern, 
                                       nChildMatern, nChild, 
                                       nGrandchildInBranches)
    # Sum array of alleles across the alleles
    alleleSumMat = apply(alleleArray, 1, function(x) {
      apply(x, 2, sum)
    })
    
    # Flag to see if any family members have pathogenic variants on both 
    # alleles for a given gene
    containsHomo = any(colSums(alleleSumMat) > 2)
    # Increment counter
    counter = counter + 1
    # Print message
    print("Simulated allele array contains pathogenic variants on both alleles for a gene, retrying.")
  }
  
  # Raise error if an admissible allele array could not be simulated 
  # after maxTries attempts
  if (counter == maxTries) {
    stop(print(paste("Could not simulate an admissible allele array for all family members after", 
                     maxTries, "attempts.")))
  }
  
  
  # Return named matrix of genotypes 
  genoMat = matrix(t(alleleSumMat != 0)*1, ncol=length(alleleFreq))
  if (is.null(names(alleleFreq))) {
    colnames(genoMat) = letters[1:ncol(genoMat)]
  } else {
    colnames(genoMat) = names(alleleFreq)
  }
  
  return(genoMat)
}