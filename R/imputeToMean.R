#' @title Impute missing values to the mean dosage
#' 
#' @description This function imputes missing genotype values to the mean
#' dosage.
#' 
#' @param Z A genotype matrix (dosage matrix) - rows correspond to 
#' individuals and columns correspond to SNPs. Use 'NA' for missing values. 
#' 
#' @details This function imputes missing genotypes to the SNP mean allele
#' frequency.
#' 
#' If Z is not a matrix it will be coerced with as.matrix. 
#' 
#' @return a genotype matrix (dosage matrix)
#'   
#' @export
imputeToMean <- function(Z) {
  
  if (!is.matrix(Z)) {
    Z <- as.matrix(Z)
  }
  
  MZ <- colMeans(Z, na.rm=TRUE)
  all0 <- which(is.nan(MZ))
  Z[ , all0] <- 0
  MZ[all0] <- 0
  ISNAZ <- is.na(Z) 
  idx <- which(ISNAZ)
  Z[idx] <- MZ[((idx-1)%/%nrow(Z))+1]
  
  return(Z)
}