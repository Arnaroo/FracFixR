# Script to create example data for FracFixR package
# This creates small datasets that will be included with the package

library(usethis)

# Set seed for reproducibility
set.seed(123)

# Create small example dataset for package examples
n_genes <- 100
n_samples <- 12

# Generate count matrix with realistic distribution
example_counts <- matrix(
  rnbinom(n_genes * n_samples, mu = 100, size = 10),
  nrow = n_genes,
  dimnames = list(
    paste0("Gene", 1:n_genes),
    paste0("Sample", 1:n_samples)
  )
)

# Add some structure to the data
# Make some genes enriched in certain fractions
# Genes 1-20: enriched in heavy polysomes
example_counts[1:20, c(3, 6, 9, 12)] <- example_counts[1:20, c(3, 6, 9, 12)] * 2

# Genes 21-40: depleted in heavy polysomes  
example_counts[21:40, c(3, 6, 9, 12)] <- round(example_counts[21:40, c(3, 6, 9, 12)] * 0.5)

# Create annotation data frame
example_annotation <- data.frame(
  Sample = colnames(example_counts),
  Condition = rep(c("Control", "Treatment"), each = 6),
  Type = rep(c("Total", "Light_Polysome", "Heavy_Polysome"), 4),
  Replicate = c(rep("Rep1", 3), rep("Rep2", 3), rep("Rep1", 3), rep("Rep2", 3)),
  stringsAsFactors = FALSE
)

# Create a polysome profiling example
polysome_annotation <- data.frame(
  Sample = paste0("Sample", 1:12),
  Condition = rep(c("Control", "Stress"), each = 6),
  Type = rep(c("Total", "Monosome", "Polysome"), 4),
  Replicate = c(rep("Rep1", 3), rep("Rep2", 3), rep("Rep1", 3), rep("Rep2", 3)),
  stringsAsFactors = FALSE
)

# Create a subcellular fractionation example
subcellular_annotation <- data.frame(
  Sample = paste0("Sample", 1:12),
  Condition = rep(c("WT", "Mutant"), each = 6),
  Type = rep(c("Total", "Nuclear", "Cytoplasmic"), 4),
  Replicate = c(rep("Rep1", 3), rep("Rep2", 3), rep("Rep1", 3), rep("Rep2", 3)),
  stringsAsFactors = FALSE
)

# Save the datasets
usethis::use_data(example_counts, overwrite = TRUE)
usethis::use_data(example_annotation, overwrite = TRUE)
usethis::use_data(polysome_annotation, overwrite = TRUE)
usethis::use_data(subcellular_annotation, overwrite = TRUE)

# Document the datasets
# Create R/data.R file with documentation
data_doc <- '
#\' Example RNA-seq count matrix
#\'
#\' A matrix containing simulated RNA-seq counts for 100 genes across 12 samples.
#\' The data simulates a polysome profiling experiment with two conditions
#\' (Control and Treatment) and three fractions (Total, Light_Polysome, Heavy_Polysome).
#\'
#\' @format A numeric matrix with 100 rows (genes) and 12 columns (samples):
#\' \\describe{
#\'   \\item{rows}{Gene identifiers (Gene1 to Gene100)}
#\'   \\item{columns}{Sample identifiers (Sample1 to Sample12)}
#\' }
#\' @source Simulated data generated for package examples
#\' @examples
#\' data(example_counts)
#\' dim(example_counts)
#\' head(example_counts[, 1:6])
"example_counts"

#\' Example annotation data frame
#\'
#\' A data frame containing sample annotations for the example_counts matrix.
#\' Describes the experimental design with conditions, fraction types, and replicates.
#\'
#\' @format A data frame with 12 rows and 4 columns:
#\' \\describe{
#\'   \\item{Sample}{Sample identifier matching column names in example_counts}
#\'   \\item{Condition}{Experimental condition (Control or Treatment)}
#\'   \\item{Type}{Fraction type (Total, Light_Polysome, or Heavy_Polysome)}
#\'   \\item{Replicate}{Replicate identifier (Rep1 or Rep2)}
#\' }
#\' @source Simulated data generated for package examples
#\' @examples
#\' data(example_annotation)
#\' head(example_annotation)
#\' table(example_annotation$Condition, example_annotation$Type)
"example_annotation"

#\' Polysome profiling annotation example
#\'
#\' An alternative annotation data frame for polysome profiling experiments
#\' with monosome and polysome fractions.
#\'
#\' @format A data frame with 12 rows and 4 columns:
#\' \\describe{
#\'   \\item{Sample}{Sample identifier}
#\'   \\item{Condition}{Experimental condition (Control or Stress)}
#\'   \\item{Type}{Fraction type (Total, Monosome, or Polysome)}
#\'   \\item{Replicate}{Replicate identifier (Rep1 or Rep2)}
#\' }
#\' @source Simulated data generated for package examples
#\' @examples
#\' data(polysome_annotation)
#\' head(polysome_annotation)
"polysome_annotation"

#\' Subcellular fractionation annotation example
#\'
#\' An annotation data frame for subcellular fractionation experiments
#\' with nuclear and cytoplasmic fractions.
#\'
#\' @format A data frame with 12 rows and 4 columns:
#\' \\describe{
#\'   \\item{Sample}{Sample identifier}
#\'   \\item{Condition}{Experimental condition (WT or Mutant)}
#\'   \\item{Type}{Fraction type (Total, Nuclear, or Cytoplasmic)}
#\'   \\item{Replicate}{Replicate identifier (Rep1 or Rep2)}
#\' }
#\' @source Simulated data generated for package examples
#\' @examples
#\' data(subcellular_annotation)
#\' head(subcellular_annotation)
"subcellular_annotation"
'

# Write documentation file
writeLines(data_doc, "R/data.R")

# Create a more complex dataset for testing
set.seed(456)
test_counts <- matrix(
  rnbinom(500 * 24, mu = 200, size = 5),
  nrow = 500,
  dimnames = list(
    paste0("Gene", 1:500),
    paste0("Sample", 1:24)
  )
)

test_annotation <- data.frame(
  Sample = colnames(test_counts),
  Condition = rep(c("Control", "Treatment"), each = 12),
  Type = rep(c("Total", "Light", "Medium", "Heavy"), 6),
  Replicate = rep(c("Rep1", "Rep2", "Rep3"), each = 4, times = 2),
  stringsAsFactors = FALSE
)

# Save test data (not exported with package, just for development)
save(test_counts, test_annotation, file = "tests/testthat/test_data.RData")

cat("Example datasets created successfully!\n")
cat("Files created:\n")
cat("  - data/example_counts.rda\n")
cat("  - data/example_annotation.rda\n")
cat("  - data/polysome_annotation.rda\n")
cat("  - data/subcellular_annotation.rda\n")
cat("  - R/data.R (documentation)\n")
cat("  - tests/testthat/test_data.RData (for testing)\n")