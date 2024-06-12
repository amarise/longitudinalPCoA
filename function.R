# R/create_plot.R

#' Create PCA Plot
#'
#' This function creates a PCA plot of the standardized residuals.
#'
#' @param data Data frame containing the microbiome data. Must have columns subjectid, time, and all OTU count data columns must begin with "OTU".
#' @param retain_prop Proportion of variability explained by top l kernel PCs.
#' @param covariates A character vector of covariates.
#' @param breaks E single value indicating the number of time intervals to cut or a numeric vector with break points.
#' @return A list containing a ggplot object, the full data frame with residuals, and the PCA proportions.
#' @export
create_plot <- function(data, retain_prop = 0.9, covariates = NULL, breaks) {
  library(lme4)
  library(ggplot2)
  library(tidyverse)
  library(Tjazi)

  # Extract columns that contain "OTU"
  Y <- as.matrix(data %>% select(contains("OTU")))

  # Embed microbiome data into Aitchison kernel matrix and get kernel PCs
  K <- AitchisonKernel(Y)

  # Retain top l kernel PCs
  l <- which.min(K$PCA_prop <= retain_prop)
  K_pcs <- K$PCs[, 1:l]
  colnames(K_pcs) <- paste0('kPC', 1:l)

  # Create a dummy data frame for example purposes
  # Replace this with your actual data frame
  # Assume 'subject', 'time', and 'response' columns in your data frame
  df <- data %>%
    select(-contains("OTU")) %>%
    cbind(K_pcs) %>%
    mutate(subjectid = factor(subjectid))

  covariate_formula <- if (is.null(covariates) || length(covariates) == 0) {
    ""
  } else {
    paste(covariates, collapse = " + ")
  }

  # Regress each kernel PC on the covariates and get estimated standardized residuals
  residuals <- sapply(1:l, function(i) {
    fit <- lmer(as.formula(paste0("kPC", i, " ~ ", covariate_formula, " + (1 | subjectid)")), data = df)
    if (isSingular(fit))
    {
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

  # Get PCs of the estimated standardized residual matrix
  residual_pca <- prcomp(residuals)

  df <- df %>%
    select(-contains('kPC')) %>%
    mutate(interval = cut(df$time, breaks = breaks, labels = FALSE)) %>%
    cbind(residual_pca$x)

  p <- ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(size = 2) +
    facet_wrap(~factor(interval), scales = 'free') +
    labs(title = "PC2 vs PC1 of Standardized Residuals",
         x = "PC1",
         y = "PC2") +
    theme_minimal()
  base_plot <- ggplot(df, aes(x = PC1, y = PC2)) +
    facet_wrap(~factor(interval), scales = 'free') +
    labs(title = "PC2 vs PC1 of Standardized Residuals",
         x = "PC1",
         y = "PC2") +
    theme_minimal()

  return(list(plot = p,
              base_plot = base_plot,
              full_df = df,
              PCA_prop = residual_pca$sdev^2/sum(residual_pca$sdev^2)))
}

#' Compute Aitchison Kernel
#'
#' This function computes the Aitchison kernel matrix from microbiome data.
#'
#' @param Y Matrix of microbiome data with samples across rows and OTUs across columns.
#' @return A list containing the kernel matrix, PCs, and PCA proportions.
#' @export
AitchisonKernel <- function(Y) {
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
