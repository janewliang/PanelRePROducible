#' Convert a Person's Allele Matrix to a Vector of Genotypes
#' 
#' Uses Mendelian laws to convert a person's alleles for a set of genes into 
#' their corresponding genotypes. Returns a vector of zeros and ones of the 
#' same length as \code{number of genes}. A 0 genotype corresponds to two 0 
#' alleles  and a 1 genotype corresponds to at least one 1 allele. 
#' @param alleles allele matrix (\code{2} by \code{number of genes})
#' @family simulations
sim.allele2Geno = function(alleles) {
  # If both alleles are 0, genotype is 0
  # Otherwise, genotype is 1 
  return(as.numeric(apply(alleles, 2, sum)!=0))
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
#' @family simulations
#' @importFrom abind abind
sim.buildGenoMat = function(alleleFreq, nChildPatern, nChildMatern, 
                            nChild, nGrandchildInBranches){ 
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
  alleles = abind(allelesPatern, allelesMatern, allelesChild, 
                  allelesGrandchildInBranches, along=1)
  
  # Return named matrix of genotypes 
  genoMat = matrix(t(apply(alleles, 1, sim.allele2Geno)), ncol=length(alleleFreq))
  if (is.null(names(alleleFreq))) {
    colnames(genoMat) = letters[1:ncol(genoMat)]
  } else {
    colnames(genoMat) = names(alleleFreq)
  }
  return(genoMat)
}