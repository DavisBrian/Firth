#' @title Run the single snp Firth test
#' 
#' @description This function runs the single snp Firth test.
#' 
#' @param formula an object of class "formula" (or one that can be coerced to 
#' that class): a symbolic description of the model to be fitted.The details 
#' of model specification are given under 'Details' in the "lm" help file.
#' @param data a data frame, list or environment (or object coercible by 
#' as.data.frame to a data frame) containing the variables in the model.
#' @param Z A genotype matrix (dosage matrix) - rows correspond to 
#' individuals and columns correspond to SNPs. Use 'NA' for missing values. 
#' The row names of this matrix should correspond to the subject names in the 
#' phenotype file.  The column names of this matrix should correspond to SNP 
#' names in the SNP information file. 
#' @param snpinfo (optional) SNP Info data.frame.  See Details.
#' @param snpNames (optional) The column in snpinfo where the SNP identifiers 
#' are found. See Details.
#' @param refCol  (optional) The column in snpinfo where the reference allele 
#' for the SNP identifier is found.  See Details.
#' @param altCol  (optional) The column in snpinfo where the alternate allele 
#' for the SNP identifier is found.  See Details.
#' 
#' @details This function is a wrapper for the logistf function in the logistf
#' package used in the Firth test. \code{\link[logistf]{logistf}}
#' 
#' @return Results of the Firth test in a RAREMETAL friendly format.
#'   
#' @export
singlesnpFirth <- function(formula, data, Z, snpinfo=NULL, snpNames=NULL, refCol=NULL, altCol=NULL) {
  
  if (!is.null(snpinfo)) {
    if (!is.null(snpNames)) {
      if (!(snpNames %in% colnames(snpinfo))) {
        stop("snpNames not a column in snpinfo")
      }
    }
    if (!is.null(refCol)) {
      if (!(refCol %in% colnames(snpinfo))) {
        stop("refCol not a column in snpinfo")
      }
    }
    if (!is.null(altCol)) {
      if (!(altCol %in% colnames(snpinfo))) {
        stop("altCol not a column in snpinfo")
      }
    }
  }
  
  if (!is.matrix(Z)) {
    Z <- as.matrix(Z)
  }
  
  mf <- model.frame(formula, data)
  
  subjects <- intersect(rownames(mf), rownames(Z))
  if (length(subjects) == 0L) {
    stop("Subject names in the phenotype and snp data do not match!")
  }
  
  Z <- imputeToMean(Z[subjects, ])
  
  y <- model.response(mf)
  cases <- names(y)[y == 1]
  controls <- names(y)[y == 0] 
  
  cols <- c("snpName", "p", "beta", "se", "MAF", "MAC", "MACcase", "MACcontrol", ".firthName")
  res <- data.frame(matrix(NA, nrow=ncol(Z), ncol=length(cols), dimnames=list(colnames(Z), cols)))
  res[ , "snpName"] <- colnames(Z)
  res[ , ".firthName"] <- paste0("firthName", 1:length(colnames(Z)))
  res[ , "MAF"] <- colMeans(Z, na.rm=TRUE)/2
  res[ , "MAC"] <- colSums(Z, na.rm=TRUE)
  res[ , "MACcase"] <- colSums(Z[cases, ], na.rm=TRUE)
  res[ , "MACcontrol"] <- colSums(Z[controls, ], na.rm=TRUE) 
  
  colnames(Z) <-  res[ , ".firthName"]
  
  mac_idx <- which(res$MAC > 0L)
  
  dat <- cbind(mf, Z[ , mac_idx])
  snp_idx <- ncol(mf)+1
  
  for (i in mac_idx) {
    f1 <- logistf::logistf(update.formula(formula, paste("~ .", colnames(Z)[i], sep="+")), data=dat)
    res[i, "beta"] <- f1$coefficients[snp_idx]
    res[i, "se"] <- sqrt(f1$var[snp_idx, snp_idx])
    #    res[i, "p"] <- f1$prob[snp_idx]
  }
  
  res$p <- 2 * pnorm(-abs(res$beta/res$se))
  
  if (!is.null(snpinfo)) {
    if (!is.null(snpNames)) {
      res <- merge(res, snpinfo[ , c(snpNames, refCol, altCol), drop=FALSE], by.x="snpName", by.y=snpNames, all=TRUE, sort=FALSE)
    }
  }
  
  return(res) 
}