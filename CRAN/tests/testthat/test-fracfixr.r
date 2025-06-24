# Test file for FracFixR package
# Run tests with: devtools::test()

context("FracFixR main functions")

# Helper function to create test data
create_test_data <- function(n_genes = 100, n_samples = 12, seed = 123) {
  set.seed(seed)
  
  # Create count matrix
  counts <- matrix(
    rnbinom(n_genes * n_samples, mu = 100, size = 10),
    nrow = n_genes,
    dimnames = list(
      paste0("Gene", 1:n_genes),
      paste0("Sample", 1:n_samples)
    )
  )
  
  # Create annotation
  annotation <- data.frame(
    Sample = paste0("Sample", 1:12),
    Condition = rep(c("Control", "Treatment"), each = 6),
    Type = rep(c("Total", "Light", "Heavy"), 4),
    Replicate = c(rep("Rep1", 3), rep("Rep2", 3), rep("Rep1", 3), rep("Rep2", 3)),
    stringsAsFactors = FALSE
  )
  
  list(counts = counts, annotation = annotation)
}

# Test FracFixR main function
test_that("FracFixR handles basic input correctly", {
  # Create test data
  test_data <- create_test_data()
  
  # Run FracFixR
  expect_silent(
    result <- FracFixR(test_data$counts, test_data$annotation)
  )
  
  # Check result structure
  expect_is(result, "list")
  expect_equal(
    names(result), 
    c("OriginalData", "Annotation", "Propestimates", "NewData", 
      "Coefficients", "Fractions", "plots")
  )
  
  # Check data dimensions
  expect_equal(nrow(result$OriginalData), nrow(test_data$counts))
  expect_equal(ncol(result$NewData), ncol(test_data$counts))
  
  # Check that fractions sum to <= 1
  frac_sums <- rowSums(result$Fractions[, !names(result$Fractions) %in% "Replicate"])
  expect_true(all(frac_sums <= 1.01)) # Allow small numerical error
})

test_that("FracFixR requires Total samples", {
  test_data <- create_test_data()
  
  # Remove Total samples
  bad_annotation <- test_data$annotation
  bad_annotation$Type[bad_annotation$Type == "Total"] <- "Other"
  
  expect_error(
    FracFixR(test_data$counts, bad_annotation),
    "at least one sample with Type == 'Total'"
  )
})

test_that("FracFixR validates input structure", {
  test_data <- create_test_data()
  
  # Test with data.frame instead of matrix
  expect_error(
    FracFixR(as.data.frame(test_data$counts), test_data$annotation),
    "must be a numeric matrix"
  )
  
  # Test with missing column names
  bad_counts <- test_data$counts
  colnames(bad_counts) <- NULL
  expect_error(
    FracFixR(bad_counts, test_data$annotation),
    "must have column names"
  )
  
  # Test with missing required annotation columns
  bad_annotation <- test_data$annotation[, -which(names(test_data$annotation) == "Type")]
  expect_error(
    FracFixR(test_data$counts, bad_annotation),
    "must contain required columns"
  )
})

test_that("FracFixR handles edge cases", {
  # Test with minimal data (2 samples per condition)
  test_data <- create_test_data(n_genes = 50, n_samples = 6)
  test_data$annotation <- data.frame(
    Sample = paste0("Sample", 1:6),
    Condition = rep(c("Control", "Treatment"), each = 3),
    Type = rep(c("Total", "Light", "Heavy"), 2),
    Replicate = rep("Rep1", 6),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    result <- FracFixR(test_data$counts, test_data$annotation),
    "fewer than 2 replicates"
  )
})

# Test ProcessReplicate function
test_that("ProcessReplicate works correctly", {
  # Create simple test data
  set.seed(123)
  test_mat <- data.frame(
    Total = c(100, 200, 150, 180, 120),
    Fraction1 = c(30, 60, 45, 54, 36),
    Fraction2 = c(50, 100, 75, 90, 60)
  )
  rownames(test_mat) <- paste0("Gene", 1:5)
  
  transcript_list <- paste0("Gene", 1:5)
  
  # Run ProcessReplicate
  result <- ProcessReplicate(test_mat, transcript_list)
  
  # Check result structure
  expect_is(result, "list")
  expect_equal(
    names(result),
    c("Propestimates", "NewData", "Coefficients", "Fractions", "plot")
  )
  
  # Check dimensions
  expect_equal(dim(result$Propestimates), dim(test_mat))
  expect_equal(dim(result$NewData), dim(test_mat))
})

test_that("ProcessReplicate requires Total column", {
  # Test data without Total column
  bad_mat <- data.frame(
    Fraction1 = c(30, 60, 45),
    Fraction2 = c(50, 100, 75)
  )
  
  expect_error(
    ProcessReplicate(bad_mat, rownames(bad_mat)),
    "must include a column named 'Total'"
  )
})

# Test DiffPropTest function
test_that("DiffPropTest performs differential testing", {
  # Create and process test data
  test_data <- create_test_data(n_genes = 50)
  fracfixr_result <- FracFixR(test_data$counts, test_data$annotation)
  
  # Test GLM method
  diff_result <- DiffPropTest(
    fracfixr_result,
    Conditions = c("Control", "Treatment"),
    Types = "Heavy",
    Test = "GLM"
  )
  
  # Check result structure
  expect_is(diff_result, "data.frame")
  expect_equal(
    names(diff_result),
    c("transcript", "mean_success_cond1", "mean_success_cond2", 
      "mean_diff", "log2FC", "pval", "padj")
  )
  
  # Check that p-values are in valid range
  valid_pvals <- !is.na(diff_result$pval)
  expect_true(all(diff_result$pval[valid_pvals] >= 0))
  expect_true(all(diff_result$pval[valid_pvals] <= 1))
  
  # Check that adjusted p-values are >= raw p-values
  valid_both <- !is.na(diff_result$pval) & !is.na(diff_result$padj)
  expect_true(all(diff_result$padj[valid_both] >= diff_result$pval[valid_both]))
})

test_that("DiffPropTest handles multiple fraction types", {
  test_data <- create_test_data(n_genes = 50)
  fracfixr_result <- FracFixR(test_data$counts, test_data$annotation)
  
  # Test with combined fractions
  diff_result <- DiffPropTest(
    fracfixr_result,
    Conditions = c("Control", "Treatment"),
    Types = c("Light", "Heavy"),
    Test = "GLM"
  )
  
  expect_is(diff_result, "data.frame")
  expect_equal(nrow(diff_result), nrow(fracfixr_result$OriginalData))
})

test_that("DiffPropTest validates inputs", {
  test_data <- create_test_data()
  fracfixr_result <- FracFixR(test_data$counts, test_data$annotation)
  
  # Test with invalid condition
  expect_error(
    DiffPropTest(fracfixr_result, 
                Conditions = c("Control", "Invalid"),
                Types = "Heavy",
                Test = "GLM"),
    "not found in NormObject\\$Annotation"
  )
  
  # Test with invalid type
  expect_error(
    DiffPropTest(fracfixr_result,
                Conditions = c("Control", "Treatment"),
                Types = "Invalid",
                Test = "GLM"),
    "not found in NormObject\\$Annotation"
  )
})

# Test visualization functions
test_that("PlotFractions creates valid plot", {
  test_data <- create_test_data()
  fracfixr_result <- FracFixR(test_data$counts, test_data$annotation)
  
  plot_obj <- PlotFractions(fracfixr_result)
  
  expect_is(plot_obj, "ggplot")
  expect_is(plot_obj, "gg")
})

test_that("PlotComparison creates valid volcano plot", {
  test_data <- create_test_data(n_genes = 50)
  fracfixr_result <- FracFixR(test_data$counts, test_data$annotation)
  
  diff_result <- DiffPropTest(
    fracfixr_result,
    Conditions = c("Control", "Treatment"),
    Types = "Heavy",
    Test = "GLM"
  )
  
  plot_obj <- PlotComparison(
    diff_result,
    Conditions = c("Control", "Treatment"),
    Types = "Heavy"
  )
  
  expect_is(plot_obj, "ggplot")
})

# Test reproducibility
test_that("FracFixR is reproducible with same seed", {
  test_data1 <- create_test_data(seed = 42)
  test_data2 <- create_test_data(seed = 42)
  
  result1 <- FracFixR(test_data1$counts, test_data1$annotation)
  result2 <- FracFixR(test_data2$counts, test_data2$annotation)
  
  # Check that coefficients are identical
  expect_equal(result1$Coefficients, result2$Coefficients)
  expect_equal(result1$Fractions, result2$Fractions)
})

# Test performance with larger datasets
test_that("FracFixR handles larger datasets", {
  skip_on_cran() # Skip on CRAN to save time
  
  # Create larger dataset
  large_data <- create_test_data(n_genes = 1000, n_samples = 24)
  
  # Time the execution
  start_time <- Sys.time()
  result <- FracFixR(large_data$counts, large_data$annotation)
  end_time <- Sys.time()
  
  # Check it completes in reasonable time (< 60 seconds)
  time_diff <- as.numeric(difftime(end_time, start_time, units = "secs"))
  expect_lt(time_diff, 60)
  
  # Check result is valid
  expect_is(result, "list")
  expect_equal(nrow(result$OriginalData), 1000)
})
