#' Example RNA-seq count matrix
#'
#' A matrix containing simulated RNA-seq counts for 100 genes across 12 samples.
#' The data simulates a polysome profiling experiment with two conditions
#' (Control and Treatment) and three fractions (Total, Light_Polysome, Heavy_Polysome).
#'
#' @format A numeric matrix with 100 rows (genes) and 12 columns (samples):
#' \describe{
#'   \item{rows}{Gene identifiers (Gene1 to Gene100)}
#'   \item{columns}{Sample identifiers (Sample1 to Sample12)}
#' }
#' @source Simulated data generated for package examples
#' @examples
#' data(example_counts)
#' dim(example_counts)
#' head(example_counts[, 1:6])
"example_counts"

#' Example annotation data frame
#'
#' A data frame containing sample annotations for the example_counts matrix.
#' Describes the experimental design with conditions, fraction types, and replicates.
#'
#' @format A data frame with 12 rows and 4 columns:
#' \describe{
#'   \item{Sample}{Sample identifier matching column names in example_counts}
#'   \item{Condition}{Experimental condition (Control or Treatment)}
#'   \item{Type}{Fraction type (Total, Light_Polysome, or Heavy_Polysome)}
#'   \item{Replicate}{Replicate identifier (Rep1 or Rep2)}
#' }
#' @source Simulated data generated for package examples
#' @examples
#' data(example_annotation)
#' head(example_annotation)
#' table(example_annotation$Condition, example_annotation$Type)
"example_annotation"

#' Polysome profiling annotation example
#'
#' An alternative annotation data frame for polysome profiling experiments
#' with monosome and polysome fractions.
#'
#' @format A data frame with 12 rows and 4 columns:
#' \describe{
#'   \item{Sample}{Sample identifier}
#'   \item{Condition}{Experimental condition (Control or Stress)}
#'   \item{Type}{Fraction type (Total, Monosome, or Polysome)}
#'   \item{Replicate}{Replicate identifier (Rep1 or Rep2)}
#' }
#' @source Simulated data generated for package examples
#' @examples
#' data(polysome_annotation)
#' head(polysome_annotation)
"polysome_annotation"

#' Subcellular fractionation annotation example
#'
#' An annotation data frame for subcellular fractionation experiments
#' with nuclear and cytoplasmic fractions.
#'
#' @format A data frame with 12 rows and 4 columns:
#' \describe{
#'   \item{Sample}{Sample identifier}
#'   \item{Condition}{Experimental condition (WT or Mutant)}
#'   \item{Type}{Fraction type (Total, Nuclear, or Cytoplasmic)}
#'   \item{Replicate}{Replicate identifier (Rep1 or Rep2)}
#' }
#' @source Simulated data generated for package examples
#' @examples
#' data(subcellular_annotation)
#' head(subcellular_annotation)
"subcellular_annotation"