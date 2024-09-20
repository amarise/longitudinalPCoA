# longitudinalPCoA

`longitudinalPCoA` is an R package designed for visualizing longitudinal and hierarchical microbiome data. It integrates linear mixed models (LMMs) and kernel matrices to accurately visualize temporal patterns, while adjusting for repeated measures and covariates.

## Features

- **Handles longitudinal data:** Incorporates linear mixed models (LMMs) to account for repeated measures.

- **Flexible visualization:** Allows customization of principal coordinate analysis (PCoA) plots.
  
- **Kernel-based approach:** Uses Aitchison (or other) kernel matrices to capture global microbiome data properties.

- **Example dataset included:** Comes with a simulated dataset for testing and demonstration.

## Installation

You can install the package directly from GitHub using devtools:

``` {r}
# Install the devtools package if necessary
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install longitudinalPCoA from GitHub
devtools::install_github("amarise/longitudinalPCoA")
```

## Getting Started

### Loading the Package and Example Data

Once installed, load the package and access the example dataset:

```{r}
# Load the package
library(longitudinalPCoA)

# Load the example dataset
data(example_data)

# View the first few rows of the metadata
head(example_data$metadata)

# Check the OTU counts matrix dimensions
dim(example_data$otu_counts)
```

The example dataset (`example_data`) contains:

- `metadata`: Information about the samples, including patient ID, time points, batch ID, and treatment group.

- `otu_counts`: A matrix of OTU counts where rows represent samples and columns represent different OTUs.

## Running Principal Coordinate Analysis (PCoA) for Repeated Measures

Here’s an example of how to perform a PCoA on the OTU counts using the Aitchison kernel matrix:

```{r}
# Perform PCoA on the OTU counts
Aitchison_kernel_matrix <- AitchisonKernel(example_data$otu_counts)$K
pcoa_result <- create_plot(
  metadata = example_data$metadata, 
  kernel_matrix = Aitchison_kernel_matrix, 
  retain_prop = 0.9, 
  covariates = c("batch"), 
  breaks = 4
)

# View the PCoA plot
pcoa_result$plot
```

## Customizing the PCoA Plot

You can customize the resulting PCoA plot by coloring the points according to treatment groups, for example:

```{r}
pcoa_result$full_df %>%
  mutate(Arm = ifelse(treatment == 0, "Control", "Treatment")) %>%
  ggplot(aes(x = PC1, y = PC2, color = Arm)) +
  geom_point(size = 2) +
  facet_wrap(~factor(interval), scales = 'free') +
  labs(title = "PC2 vs PC1 of Standardized Residuals",
         x = "PC1",
         y = "PC2") +
  theme_minimal()
```

Documentation and Vignettes

For a detailed guide on using the package, see the vignette:

```{r}
browseVignettes("longitudinalPCoA")
```

The vignette provides step-by-step instructions on loading data, running analyses, and customizing visualizations.

## Conclusion

`longitudinalPCoA` provides a robust framework for visualizing microbiome data over time. By accounting for repeated measures and allowing flexible customization of PCoA plots, the package helps researchers gain clearer insights into temporal microbiome dynamics.
