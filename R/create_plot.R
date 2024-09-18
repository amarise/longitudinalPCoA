#' Create PCA Plot
#'
#' This function creates a PCA plot of the standardized residuals.
#'
#' @param data Data frame containing the microbiome data. Must have columns subjectid, time, and all OTU count data columns must begin with "OTU".
#' @param retain_prop Proportion of variability explained by top l kernel PCs.
#' @param covariates A character vector of covariates.
#' @param breaks A single value indicating the number of time intervals to cut or a numeric vector with break points.
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
      var.y <- var.b + vr * Diagonal(nrow(df))
      cholVarInv <- solve(t(chol(var.y)))
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
