# usethis::create_package("path/to/longitudinalPCoA")
usethis::use_readme_md()
# usethis::use_data(example_data)
usethis::use_mit_license()
usethis::use_testthat()
usethis::use_vignette("Introduction-to-longitudinalPCoA")
usethis::use_build_ignore(c("README.md", "tests"))
usethis::use_git_ignore(c(".Rhistory", ".RData", ".Rproj.user", ".DS_Store", ".Rbuildignore"))
library(testthat)
library(longitudinalPCoA)

load("~/Documents/Longitudinal Microbiome Visualization/example_data.RDa")

# Save the dataset for use in the package
usethis::use_data(example_data, overwrite = TRUE)
