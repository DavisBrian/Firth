#' @title Code to the Minor Allele
#' 
#' @description This function ensures the genotype is coded to the minor allele.
#' 
#' @param Z A genotype matrix (dosage matrix) - rows correspond to 
#' individuals and columns correspond to SNPs. Use 'NA' for missing values. 
#' 
#' @details This function ensures the minor allele frequency is < 0.5.  Note 
#' this assumes the genotype matrix is coded 0/1/2.
#' 
#' If Z is not a matrix it will be coerced with as.matrix. 
#' 
#' @return a recoded genotype matrix (dosage matrix)
#'   
#' @export
codeToMinor <- function(Z) {
  
  if (!is.matrix(Z)) {
    Z <- as.matrix(Z)
  }
  
  maf <- colMeans(Z, na.rm=TRUE)/2.0
  idx <- which(maf > .50)
  Z[ , idx] <- 2 - Z[ , idx]
  
  return(Z)
}
