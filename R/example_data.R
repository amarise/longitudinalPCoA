#' Example Dataset: Simulated Microbiome Data
#'
#' This dataset is a list containing simulated microbiome data from a study with 100 patients, each measured at 4 time points. It includes sample metadata and OTU (Operational Taxonomic Unit) counts.
#'
#' @format A list with 2 elements:
#' \describe{
#'   \item{metadata}{A data frame containing sample metadata with 400 rows (100 patients, each measured 4 times). Columns include:
#'     \itemize{
#'       \item \code{subjectid}: The ID of each subject
#'       \item \code{time}: The time point of sample collection.
#'       \item \code{batch}: The batch ID indicating different experimental batches.
#'       \item \code{treatment}: Binary treatment indicator (0 or 1) indicating the treatment group of the patient.
#'     }
#'   }
#'   \item{otu_counts}{A matrix of OTU counts with 400 rows (samples) and 233 columns (OTUs). Each value corresponds to the observed count of a specific OTU in a given sample.}
#' }
#'
#' @examples
#' data(example_data)
#' str(example_data)  # View the structure of the dataset
#' example_data$metadata  # Access the metadata data frame
#' example_data$otu_counts  # Access the OTU counts matrix with 233 OTUs
"example_data"
