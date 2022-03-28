#' Simulate Current Ages for Branch
#'
#' Returns a numeric vector of current ages for the progeny of the  
#' branch, plus any parents who have missing ages (parents who were simulated/not 
#' blood relations of the main family, i.e. spouses or the grandparents). The 
#' childrens' ages are simulated based on their mother's age. 
#' @param nProg number of progeny in branch
#' @param mothAge mother's age, if known (optional)
#' @param fathAge father's age, if known (optional)
#' @details Does not ensure that the children have different ages from each other
#' (i.e., two siblings with the same ages are not necessarily twins).
#' If one of the parents' ages is missing, it will be simulated based on the other 
#' parent's age. If both are missing, they will be assumed to be at the top of tree 
#' (grandparents). The grandmother's age will be simulated from a Norm(100, 2) 
#' distribution and the grandfather's age will be simulated based on her age. 
#' @family simulations
sim.simCurAgeVar = function(nProg, mothAge, fathAge) {
  # If both parents' ages are missing, assume they are at the top of the tree 
  # (grandparents)
  if (missing(mothAge) & missing(fathAge)) {
    # Simulate both parents' ages and add them to the vector of current ages
    # Center mother's age distribution at 100 
    # Center father's age distribution around mother's age
    mothAge = round(rnorm(1, 100, 2))
    fathAge = round(rnorm(1, mothAge, 2))
    CurAge = c(fathAge, mothAge)
    # If only the mother's age is missing, simulate it to be close to the father's age
    # Add it to the vector of current ages
    # Must be at least 15
  } else if (missing(mothAge)) {
    mothAge = max(15, round(rnorm(1, fathAge, 2)))
    CurAge = mothAge
    # If only the father's age is missing, simulate it to be close to the mother's age
    # Add it to the vector of current ages
    # Must be at least 15
  } else if (missing(fathAge)) {
    fathAge = max(15, round(rnorm(1, mothAge, 2)))
    CurAge = fathAge
  } else {
    CurAge = numeric(0)
  }
  
  # Simulate the children's ages to be centered around the mother's age minus 30
  CurAge = c(CurAge, pmin(mothAge-15, round(rnorm(nProg, mothAge-30, 5))))
  return(CurAge)
}