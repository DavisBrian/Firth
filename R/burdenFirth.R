#' @title Run the burden Firth test
#' 
#' @description This function runs the burden Firth test.
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
#' @param snpinfo SNP Info data.frame. Must contain fields given in 'snpName' 
#' and 'aggregateBy'
#' @param snpNames The column in snpinfo where the SNP identifiers.  
#' Defualt is 'Name'
#' @param aggregateBy  The column in snpinfo on which the individual snps are 
#' to be aggregated by. Default is 'gene'. 
#' @param mafRange  Range of MAF's to include in the analysis (endpoints 
#' included). Default is all SNPs (0 <= MAF <= 0.5).
#' 
#' @details This function is a wrapper for the logistf function in the logistf
#' package used in the Firth test.
#' 
#' @return Results of the burden Firth test.
#'   
#' @export
burdenFirth <- function(formula, data, Z, snpinfo, aggregateBy="gene", snpNames="Name", mafRange = c(0,0.5)) {
  
  if (!is.matrix(Z)) {
    Z <- as.matrix(Z)
  }
  
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  cases <- names(y)[y == 1]
  controls <- names(y)[y == 0] 
  
  subjects <- intersect(names(y), rownames(Z))
  if (length(subjects) == 0L) {
    stop("Subject names in the phenotype and snp data do not match!")
  }
  
  snps <- intersect(colnames(Z), snpinfo[ , snpNames])
  if (length(snps) == 0L) {
    stop("SNP names in the snpinfo and snp data do not match!")
  }
  
  #  Z <- Z[subjects, snps]
  Z <- imputeToMean(codeToMinor(Z[subjects, snps]))
  
  snp.name.list <- split(snpinfo[ , snpNames], snpinfo[ , aggregateBy])
  genes <- names(snp.name.list)
  
  cols <- c("Gene", "p", "beta", "se", "NSNPS", "NSNPSused", "MAC", "MACcase", "MACcontrol", "doseMean", "doseVar")
  res <- data.frame(matrix(NA, nrow=length(genes), ncol=length(cols), dimnames=list(genes, cols)))
  
  dose <- matrix(0, nrow=nrow(Z), ncol=length(genes), dimnames=list(rownames(Z), genes))
  
  maf <- colMeans(Z, na.rm=TRUE)/2
  keep <- (maf > min(mafRange)) & (maf <= max(mafRange))
  keep[is.na(keep)] <- FALSE
  keep.snps <- names(keep)[keep]  
  
  for (gene in genes) {
    res[gene, "Gene"] <- gene
    gene.snps.all <- intersect(snp.name.list[[gene]], colnames(Z))
    gene.snps.keep <- intersect(gene.snps.all, keep.snps)
    res[gene, "NSNPS"] <- length(gene.snps.all)
    res[gene, "NSNPSused"] <- length(gene.snps.keep)
    
    if (length(gene.snps.keep) > 0L){
      res[gene, "MAC"] <- sum(Z[ , gene.snps.keep], na.rm=TRUE)
      res[gene, "MACcase"] <- sum(Z[cases , gene.snps.keep], na.rm=TRUE)
      res[gene, "MACcontrol"] <- sum(Z[controls , gene.snps.keep], na.rm=TRUE)    
      
      if (length(gene.snps.keep) == 1L) {
        dose[subjects , gene] <- Z[subjects , gene.snps.keep]      
      } else {
        dose[subjects, gene] <- rowSums(Z[subjects , gene.snps.keep], na.rm=TRUE)
      }
      # can't use colSums, colMeans, etc
      res[gene, "doseMean"] <- mean(dose[subjects, gene], na.rm=TRUE)
      res[gene, "doseVar"] <- var(dose[subjects, gene], na.rm=TRUE)
    } 
  }
  
  tmp <- singlesnpFirth(formula=formula, data=data, Z=dose)
  
  # merge
  gene.list <- intersect(rownames(res), rownames(tmp))
  res[gene.list, c("p", "beta", "se")] <- tmp[gene.list, c("p", "beta", "se")]
  
  return(res)
}
