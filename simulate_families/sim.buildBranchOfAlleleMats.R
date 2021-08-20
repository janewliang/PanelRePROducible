#' Simulate a Person's Matrix of Alleles for a Series of Genes
#' 
#' Returns a \code{2} by \code{number of genes} matrix of zeros and ones, 
#' indicating which alleles the person has for each of the \code{length(alleleFreq)} 
#' genes, simulated based on the supplied allele frequencies.
#' @param alleleFreq allele frequencies for each gene of interest 
#' @family simulations
sim.simAlleleMat = function(alleleFreq) {
  return(sapply(alleleFreq, function(x){replicate(2, rbinom(1, 1, x))}))
}


#' Simulate a Person's Matrix of Alleles for a Series of Genes, Based on Inheritance
#' 
#' Simulate a person's alleles for a series of genes, based on their parents' alleles 
#' and Mendelian laws of inheritance. Returns a \code{2} by \code{number of genes} 
#' matrix of zeros and ones, where the first row corresponds to the alleles inherited
#' from the father and the second row corresponds to the alleles inherited from the
#' mother.  
#' @param alleles1 \code{2} by \code{number of genes} allele matrix for parent 1 
#' @param alleles2 \code{2} by \code{number of genes} allele matrix for parent 2  
#' @family simulations
sim.inheritAlleleMat = function(alleles1, alleles2) {
  # Number of genes
  nGenes = ncol(alleles1)
  if (nGenes != ncol(alleles2)) {
    stop("The parents have allele matrices for a different number of genes.")
  }
  
  # Sample the allele for each gene that the person will inherit
  idx1 = sample(2, nGenes, replace=TRUE) # Parent 1's alleles
  idx2 = sample(2, nGenes, replace=TRUE) # Parent 2's alleles
  
  # Person's allele matrix
  # First row corresponds to alleles inherited from father
  # Second row corresponds to alleles inherited from mother
  return(rbind(alleles1[cbind(idx1,1:nGenes)], 
               alleles2[cbind(idx2,1:nGenes)]))
}


#' Build the Allele Matrices for the Progeny of a Branch
#' 
#' Returns a 3d array of allele matrices with dimensions for the progeny (and 
#' possibly one or more of the parents) of the branch. The first dimension of the 3d 
#' array will be equal to the number of progeny plus up to two additional parents. 
#' The other two dimensions are \code{2} by \code{number of genes}. See Details. 
#' @param alleleFreq allele frequencies for each gene of interest
#' @param nProg number of progeny in this branch
#' @param alleles \code{2} by \code{number of genes} matrix of alleles 
#' for family member who creates this branch (optional)
#' @param allelesSpouse \code{2} by \code{number of genes} matrix of alleles 
#' for family member's spouse (optional)
#' @details If \code{alleles} and/or \code{allelesSpouse} are not supplied, they will 
#' be simulated. In this case, the simulated allele matrices will be added along the 
#' first dimension of the 3d array, i.e. the first dimension will grow to 
#' \code{nProg+1} or \code{nProg+2}. 
#' @family simulations
sim.buildBranchOfAlleleMats = function(alleleFreq, nProg, alleles, allelesSpouse) {
  # Check if the allele matrix for the branch head is not supplied
  isAddHead = missing(alleles)
  # Check if the allele matrix for the spouse is not supplied
  isAddSpouse = missing(allelesSpouse)
  
  # Simulate the branch head's allele matrix if not passed in 
  if (isAddHead) {
    alleles = sim.simAlleleMat(alleleFreq)
  }
  # Simulate the spouse's allele matrix if not passed in 
  if (isAddSpouse) {
    allelesSpouse = sim.simAlleleMat(alleleFreq)
  }
  
  # Generate allele matrices for each progeny, using Mendelian laws of inheritance
  if (nProg > 0) {
    allelesProg = aperm(replicate(nProg, 
                                  sim.inheritAlleleMat(alleles, allelesSpouse)),
                        c(3,1,2))
  } else {
    allelesProg = array(dim=c(nProg, dim(alleles)))
  }
  
  # Put everything into the same 3d matrix
  # If branch head's allele matrix was not supplied, include the simulated matrix
  if (isAddHead) { 
    allelesBranch = array(dim=c(nProg+2, dim(alleles)))
    allelesBranch[1,,] = alleles
    allelesBranch[2,,] = allelesSpouse
    allelesBranch[-c(1,2),,] = allelesProg
    # If spouse's allele matrix was not supplied, include the simulated matrix
  } else if (isAddSpouse) {
    allelesBranch = array(dim=c(nProg+1, dim(alleles)))
    allelesBranch[1,,] = allelesSpouse
    allelesBranch[-1,,] = allelesProg
  } else {
    allelesBranch = allelesProg
  }
  
  return(allelesBranch)
}
