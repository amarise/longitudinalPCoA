#' Compute Bray Curtis Kernel
#'
#' This function computes the Bray Curtis kernel matrix from microbiome data.
#'
#' @param Y Matrix of microbiome count data with samples across rows and OTUs across columns.
#' @return A list containing the kernel matrix, PCs, and PCA proportions.
#' @export
BrayCurtisKernel <- function(Y) {
  library(vegan)
  otu.D <- as.matrix(vegdist(Y, method="bray"))  # B-C distance

  # convert distance to kernel matrix
  n <- nrow(otu.D)
  centerM <- diag(n) - 1/n
  K <- -0.5 * centerM %*% (otu.D * otu.D) %*% centerM
  eK <- eigen(K, symmetric = TRUE)
  K <- eK$vector %*% diag(pmax(0, eK$values)) %*% t(eK$vector)

  mK <- which(eK$values > 1e-10)
  PCs <- eK$vector[,mK] %*% diag(sqrt(eK$values[mK]))
  PCA_prop <- cumsum(eK$values[mK]/sum(eK$values[mK]))

  return(list(K = K,
              PCs = PCs,
              PCA_prop = PCA_prop))
}
