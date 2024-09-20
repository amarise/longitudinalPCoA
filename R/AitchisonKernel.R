#' Compute Aitchison Kernel
#'
#' This function computes the Aitchison kernel matrix from microbiome data.
#'
#' @param Y Matrix of microbiome count data with samples across rows and OTUs across columns.
#' @return A list containing the kernel matrix, PCs, and PCA proportions.
#' @export
AitchisonKernel <- function(Y) {
  library(Tjazi)
  clr <- clr_lite(Y, samples_are = 'rows', method = 'logunif', replicates = 1000)
  n <- nrow(Y)
  centerM <- diag(n) - 1/n
  D <- as.matrix(dist(clr))
  K <- -0.5 * centerM %*% (D * D) %*% centerM
  eK <- eigen(K, symmetric = TRUE)
  K <- eK$vector %*% tcrossprod(diag(pmax(0, eK$values)), eK$vector)
  mK <- which(eK$values > 1e-10)
  PCs <- eK$vector[,mK] %*% diag(sqrt(eK$values[mK]))
  PCA_prop <- cumsum(eK$values[mK]/sum(eK$values[mK]))
  return(list(K = K,
              PCs = PCs,
              PCA_prop = PCA_prop))
}
