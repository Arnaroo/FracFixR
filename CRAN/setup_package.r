# Setup script for creating FracFixR package structure
# Run this script to set up the complete package

# Install necessary development packages if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("usethis", quietly = TRUE)) install.packages("usethis")
if (!requireNamespace("roxygen2", quietly = TRUE)) install.packages("roxygen2")
if (!requireNamespace("testthat", quietly = TRUE)) install.packages("testthat")
if (!requireNamespace("covr", quietly = TRUE)) install.packages("covr")

library(devtools)
library(usethis)

# Set your working directory to where you want to create the package
# setwd("~/projects")

# Create package structure
# usethis::create_package("FracFixR")

# Set working directory to package root
# setwd("FracFixR")

# Configure package settings
usethis::use_description(
  fields = list(
    Title = "Compositional Statistical Framework for RNA Fractionation Analysis",
    Description = "A compositional statistical framework for absolute proportion estimation between fractions in RNA sequencing data. FracFixR addresses the fundamental challenge in fractionated RNA-seq experiments where library preparation and sequencing depth obscure the original proportions of RNA fractions. It reconstructs original fraction proportions using non-negative linear regression, estimates the 'lost' unrecoverable fraction, corrects individual transcript frequencies, and performs differential proportion testing between conditions. Supports any RNA fractionation protocol including polysome profiling, subcellular localization, and RNA-protein complex isolation.",
    `Authors@R` = c(
      person("Alice", "Cleynen", email = "alice.cleynen@umontpellier.fr", 
             role = c("aut", "cre"), comment = c(ORCID = "0000-0000-0000-0000")),
      person("Agin", "Ravindran", role = "aut"),
      person("Nikolay", "Shirokikh", email = "nikolay.shirokikh@uwa.edu.au", 
             role = "aut", comment = c(ORCID = "0000-0000-0000-0000"))
    ),
    License = "CC BY-NC-ND 4.0",
    URL = "https://github.com/Arnaroo/FracFixR",
    BugReports = "https://github.com/Arnaroo/FracFixR/issues",
    Version = "1.0.0"
  )
)

# Set up version control
usethis::use_git()
usethis::use_github()  # This will create the GitHub repo if authenticated

# Set up documentation
usethis::use_readme_md()
usethis::use_news_md()

# Set up testing infrastructure
usethis::use_testthat()

# Set up continuous integration
usethis::use_github_action_check_standard()
usethis::use_github_action("test-coverage")

# Create necessary directories
dir.create("inst", showWarnings = FALSE)
dir.create("data-raw", showWarnings = FALSE)

# Set up vignette infrastructure
usethis::use_vignette("FracFixR-intro")

# Create example data
# Place this in data-raw/create_example_data.R
example_data_script <- '
# Script to create example data for FracFixR package
set.seed(123)

# Create small example dataset
n_genes <- 100
n_samples <- 12

example_counts <- matrix(
  rnbinom(n_genes * n_samples, mu = 100, size = 10),
  nrow = n_genes,
  dimnames = list(
    paste0("Gene", 1:n_genes),
    paste0("Sample", 1:n_samples)
  )
)

example_annotation <- data.frame(
  Sample = colnames(example_counts),
  Condition = rep(c("Control", "Treatment"), each = 6),
  Type = rep(c("Total", "Light_Polysome", "Heavy_Polysome"), 4),
  Replicate = c(rep("Rep1", 3), rep("Rep2", 3), rep("Rep1", 3), rep("Rep2", 3)),
  stringsAsFactors = FALSE
)

# Save the data
usethis::use_data(example_counts, example_annotation, overwrite = TRUE)
'

# Write the data creation script
writeLines(example_data_script, "data-raw/create_example_data.R")

# Run roxygen to generate documentation
roxygen2::roxygenise()

# Run tests
devtools::test()

# Check the package
devtools::check()

# Build the package
devtools::build()

# Optional: Check with stricter CRAN settings
devtools::check(cran = TRUE)

# Check on different platforms (requires internet)
# devtools::check_win_devel()
# devtools::check_win_release()
# devtools::check_mac_release()

# Calculate test coverage
# covr::package_coverage()

# Final steps before submission:
# 1. Update cran-comments.md with check results
# 2. Run spell check: devtools::spell_check()
# 3. Check URLs: urlchecker::url_check()
# 4. Run final check: devtools::check(cran = TRUE)
# 5. Submit to CRAN: devtools::release()