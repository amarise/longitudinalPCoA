#' Create PCA Plot
#'
#' This function creates a PCA plot of the standardized residuals based on metadata and a kernel matrix.
#'
#' @param metadata Data frame containing the microbiome metadata. Must have columns `subjectid`, `time`, and other covariates.
#' @param kernel_matrix A kernel matrix derived from OTU data (e.g., Aitchison kernel).
#' @param retain_prop Proportion of variability explained by top l kernel PCs. Default is 0.9.
#' @param covariates A character vector of covariates to adjust for in the linear mixed model. Must be columns in `metadata`
#' @param breaks A single value indicating the number of time intervals to cut or a numeric vector with break points.
#' @return A list containing a ggplot object, the full data frame with residuals, and the PCoA proportions.
#' @export
create_plot <- function(metadata, kernel_matrix, retain_prop = 0.9, covariates = NULL, breaks) {
  library(lme4)
  library(ggplot2)
  library(tidyverse)

  # Get the PCs from the kernel matrix
  K <- kernel_matrix
  eigen_K <- eigen(K, symmetric = TRUE)

  # Retain the top l kernel PCs
  total_variance <- cumsum(eigen_K$values / sum(eigen_K$values))
  l <- which.min(total_variance <= retain_prop)
  K_pcs <- eigen_K$vectors[, 1:l]
  colnames(K_pcs) <- paste0('kPC', 1:l)

  # Combine metadata and kernel PCs
  df <- cbind(metadata, K_pcs)
  df$subjectid <- factor(df$subjectid)

  # Create covariate formula
  covariate_formula <- if (is.null(covariates) || length(covariates) == 0) {
    ""
  } else {
    paste(covariates, collapse = " + ")
  }

  # Regress each kernel PC on the covariates and get estimated standardized residuals
  residuals <- sapply(1:l, function(i) {
    fit <- lmer(as.formula(paste0("kPC", i, " ~ ", covariate_formula, " + (1 | subjectid)")), data = df)
    if (isSingular(fit)) {
      stop(paste("The model for kPC", i, "is singular."))
    } else {
      # compute covariance matrix of random effects
      var.d <- crossprod(getME(fit,"Lambdat"))
      Zt <- getME(fit, "Zt")
      vr <- sigma(fit)^2
      var.b <- vr * (crossprod(Zt, var.d) %*% Zt)
      # compute the covariance matrix of the response variable
      var.y <- var.b + vr * Diagonal(nrow(df))
      # compute the inverse of the Cholesky decomposition
      cholVarInv <- solve(t(chol(var.y)))
      # compute standardized residuals
      as.numeric(cholVarInv %*% (df[,paste0('kPC',i)] - predict(fit, re.form = NA)))
    }
  })

  # Perform PCA on the residuals
  residual_pca <- prcomp(residuals)

  # Add PCA results to the data frame and add time intervals
  df <- df %>%
    mutate(interval = cut(df$time, breaks = breaks, labels = FALSE)) %>%
    cbind(residual_pca$x)

  # Plot PC1 vs PC2 of the standardized residuals
  p <- ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(size = 2) +
    facet_wrap(~factor(interval), scales = 'free') +
    labs(title = "PC2 vs PC1 of Standardized Residuals",
         x = "PC1",
         y = "PC2") +
    theme_minimal()

  return(list(plot = p,
              full_df = df,
              PCA_prop = residual_pca$sdev^2 / sum(residual_pca$sdev^2)))
}
