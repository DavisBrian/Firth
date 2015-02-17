#' Firth: A package for computating the Firth 
#'
#' The Firth package provides an user friendly way to run Firth's
#' penalized-likelihood logistic regression on genomic data.
#' 
#' @section Firth functions:
#' \itemize{
#'  \item \code{burdenFirth}: runs the burden Firth test
#'  
#'  \item \code{codeToMinor}: ensures the genotype matrix is coded to the minor allele
#'  
#'  \item \code{imputeToMean}: imputes missing genotype values to the mean dosage
#'  
#'  \item \code{singsnpFirth}: runs the single snp Firth test
#'  
#' }
#'
#' @docType package
#' @name Firth
#' @import logistf
NULL