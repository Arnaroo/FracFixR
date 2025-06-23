################################################################################
# FracFixR: Fraction Correction Framework for RNA-seq Data
################################################################################
# 
# Description:
#   A compositional statistical framework for absolute proportion estimation 
#   between fractions in RNA sequencing data. FracFixR addresses the fundamental
#   challenge in fractionated RNA-seq experiments where library preparation and
#   sequencing depth obscure the original proportions of RNA fractions.
#
# Key Features:
#   - Reconstructs original fraction proportions using non-negative linear regression
#   - Estimates the "lost" unrecoverable fraction
#   - Corrects individual transcript frequencies 
#   - Performs differential proportion testing between conditions
#   - Handles any RNA fractionation protocol (polysome profiling, subcellular, etc.)
#
# Authors: Alice Cleynen, Agin Ravindran, Nikolay Shirokikh
# License: [To be specified]
# Version: 1.0.0
################################################################################

# ==============================================================================
# SECTION 1: PACKAGE DEPENDENCIES
# ==============================================================================
# Load required packages for parallel processing, optimization, visualization, 
# and data manipulation

library(future.apply)      # Parallel processing capabilities
library(nnls)             # Non-negative least squares regression
library(ggplot2)          # Advanced plotting
library(dplyr)            # Data manipulation
library(RColorBrewer)     # Color palettes for plots
library(tidyr)            # Data tidying operations
library(matrixStats)      # Efficient matrix operations
library(EnhancedVolcano)  # Volcano plot visualization
library(aod)              # Beta-binomial regression for Wald test

# ==============================================================================
# SECTION 2: MAIN FRACFIXR FUNCTION
# ==============================================================================

#' FracFixR: Main Function for Fraction Correction
#'
#' @description
#' This is the core function that implements the FracFixR framework. It takes
#' raw count data from total and fractionated samples and reconstructs the
#' original fraction proportions through compositional modeling.
#'
#' @param MatrixCounts A numeric matrix of raw transcript/gene counts with:
#'   - Rows: transcripts/genes (must have rownames)
#'   - Columns: samples (must have colnames matching Annotation$Sample)
#' @param Annotation A data.frame with required columns:
#'   - Sample: sample identifiers matching column names in MatrixCounts
#'   - Condition: experimental condition (e.g., "Control", "Treatment")
#'   - Type: fraction type (must include at least one "Total")
#'   - Replicate: replicate identifier
#'
#' @return A list containing:
#'   - OriginalData: filtered input count matrix
#'   - Annotation: input annotation data
#'   - Propestimates: matrix of proportion estimates
#'   - NewData: matrix of corrected counts
#'   - Coefficients: data.frame of regression coefficients
#'   - Fractions: data.frame of estimated fraction proportions
#'   - plots: list of diagnostic plots
#'
#' @details
#' The function works by:
#' 1. Filtering transcripts based on presence in Total samples
#' 2. For each condition and replicate, fitting NNLS regression
#' 3. Estimating global fraction weights and individual transcript proportions
#' 4. Calculating the "lost" unrecoverable fraction
FracFixR <- function(MatrixCounts, Annotation) {
  # --------------------------------------------------------------------------
  # Input validation: Ensure matrix format and required structure
  # --------------------------------------------------------------------------
  stopifnot(is.matrix(MatrixCounts),
            is.numeric(MatrixCounts),
            !is.null(rownames(MatrixCounts)),
            !is.null(colnames(MatrixCounts)))
  
  stopifnot(is.data.frame(Annotation),
            all(c("Sample","Condition","Type","Replicate") %in% colnames(Annotation)))
  
  # Verify all samples in counts matrix are documented in annotation
  if (!all(colnames(MatrixCounts) %in% Annotation$Sample))
    stop("All MatrixCounts columns must be in Annotation$Sample")
  
  # Check that at least one Total sample is present (required for normalization)
  if (!any(Annotation$Type == "Total")) {
    stop("Annotation must include at least one sample with Type == 'Total'")
  }
  
  # Check for minimum replicates per condition
  reps_per_condition <- Annotation %>%
    group_by(Condition) %>%
    summarise(n_reps = n_distinct(Replicate), .groups = "drop")
  
  if (any(reps_per_condition$n_reps < 2)) {
    warning("Some conditions have fewer than 2 replicates. Statistical tests may have limited power.")
  }
  
  # --------------------------------------------------------------------------
  # Setup parallel processing for computational efficiency
  # --------------------------------------------------------------------------
  message("Setting up parallel processing...")
  future::plan("multisession", workers = parallel::detectCores() - 1)
  
  # --------------------------------------------------------------------------
  # Filter transcripts: Keep only those present in Total samples
  # This ensures we work with transcripts that are detectable in the whole cell
  # --------------------------------------------------------------------------
  message("Filtering transcripts based on Total samples...")
  wct_idx <- which(Annotation$Type == "Total")
  if (length(wct_idx) > 1) {
    # Multiple Total samples: keep transcripts present in at least one
    DataNorm <- MatrixCounts[rowSums(MatrixCounts[, wct_idx, drop = FALSE]) > 0, , drop = FALSE]
  } else {
    # Single Total sample: keep transcripts with non-zero counts
    DataNorm <- MatrixCounts[MatrixCounts[, wct_idx] > 0, , drop = FALSE]
  }
  message(sprintf("Retained %d transcripts present in Total samples", nrow(DataNorm)))
  
  # Combine transposed count data with annotation for easier manipulation
  Data <- cbind(t(DataNorm), Annotation[match(colnames(DataNorm), Annotation$Sample), ])
  
  # --------------------------------------------------------------------------
  # Initialize storage containers for results
  # --------------------------------------------------------------------------
  NewDataComplete <- NULL        # Corrected count matrix
  PropestimatesComplete <- NULL  # Proportion estimates
  CoefficientComplete <- list()  # Regression coefficients
  FractionsComplete <- list()    # Fraction proportions
  all_plots <- list()           # Diagnostic plots
  rep_ids <- character()        # Replicate identifiers

  # --------------------------------------------------------------------------
  # Process each experimental condition separately
  # --------------------------------------------------------------------------
  n_conditions <- length(unique(Annotation$Condition))
  cond_counter <- 0
  
  for (cond in unique(Annotation$Condition)) {
    cond_counter <- cond_counter + 1
    message(sprintf("\nProcessing condition %d/%d: %s", cond_counter, n_conditions, cond))
    
    # Extract data for current condition
    Datacond <- subset(Data, Condition == cond)
    DataDD <- DataNorm[, Annotation$Condition == cond, drop = FALSE]
    AnnoCond <- subset(Annotation, Condition == cond)
    wct_cond <- which(AnnoCond$Type == "Total")

    # Calculate total RNA abundance across Total samples for this condition
    if (length(wct_cond) > 1) {
      TotalSum <- rowSums(DataDD[, wct_cond, drop = FALSE])
    } else {
      TotalSum <- DataDD[, wct_cond]
    }
    if (!is.numeric(TotalSum))
      stop("TotalSum must be numeric")

    # ----------------------------------------------------------------------
    # Select informative transcripts for regression
    # Use 70-96% quantile range to avoid:
    # - Low abundance transcripts (noisy)
    # - Very high abundance transcripts (potential outliers)
    # ----------------------------------------------------------------------
    s1 <- quantile(TotalSum, 0.7)
    s2 <- quantile(TotalSum, 0.96)
    message(sprintf("  Selecting transcripts with 70-96%% quantiles (range: %.1f - %.1f)", s1, s2))

    # Create list of transcripts in the selected abundance range
    # Ensure at least 7 counts to avoid very low abundance noise
    transcriptlist <- rownames(DataDD)[TotalSum > max(s1,7) & TotalSum < s2]
    message(sprintf("  Selected %d transcripts for regression", length(transcriptlist)))

    # Check if we have enough transcripts
    if (length(transcriptlist) == 0) {
      stop(sprintf("No transcripts selected for regression in condition %s", cond))
    }
    
    if (length(transcriptlist) < 10) {
      warning(sprintf("Only %d transcripts selected for condition %s (minimum recommended: 10)", 
                     length(transcriptlist), cond))
    }

    # ----------------------------------------------------------------------
    # Process each replicate within the condition
    # ----------------------------------------------------------------------
    n_reps <- length(unique(Datacond$Replicate))
    rep_counter <- 0
    
    for (rep in unique(Datacond$Replicate)) {
      rep_counter <- rep_counter + 1
      message(sprintf("  Processing replicate %d/%d: %s", rep_counter, n_reps, rep))
      
      # Extract replicate data and ensure Total is first
      Datatemp <- subset(Datacond, Replicate == rep)
      Datatemp <- Datatemp[order(Datatemp$Type != "Total"), ]
      
      # Extract count measurements (exclude metadata columns)
      measurements <- Datatemp[, setdiff(names(Datatemp), c("Sample","Type","Replicate","Condition"))]
      if (nrow(measurements) < 2)
        stop("Each replicate must have at least two types (Total + one fraction)")
      
      # Replace first row with calculated TotalSum for consistency
      measurements[1, ] <- TotalSum
      sample_names <- Datatemp$Sample
      sample_names[1] <- "Total"  # Standardize Total sample name
      
      # Prepare data matrix for ProcessReplicate function
      DataProcess <- as.data.frame(t(measurements))
      colnames(DataProcess) <- sample_names

      # ------------------------------------------------------------------
      # Core processing: fit NNLS model for this replicate
      # ------------------------------------------------------------------
      results <- ProcessReplicate(DataProcess, transcriptlist)
      
      # Validate results structure
      stopifnot(is.list(results),
                all(c("NewData","Propestimates","Coefficients","Fractions", "plot") %in% names(results)))

      # ------------------------------------------------------------------
      # Assign meaningful names to results and compile
      # ------------------------------------------------------------------
      colnames(results$NewData) <- Datatemp$Sample
      colnames(results$Propestimates) <- Datatemp$Sample
      names(results$Coefficients) <- Datatemp$Type
      names(results$Coefficients)[1] <- "Lost"  # First coefficient is intercept (lost fraction)
      names(results$Fractions) <- Datatemp$Type[-1]  # Exclude Total from fractions

      # Bind results to complete matrices
      # Note: For very large datasets, consider pre-allocating matrices for efficiency
      NewDataComplete <- if (is.null(NewDataComplete)) results$NewData else cbind(NewDataComplete, results$NewData)
      PropestimatesComplete <- if (is.null(PropestimatesComplete)) results$Propestimates else cbind(PropestimatesComplete, results$Propestimates)

      # Store replicate-specific results with unique identifiers
      FractionsComplete[[paste(cond, rep, sep = "_")]] <- results$Fractions
      CoefficientComplete[[paste(cond,rep,sep="_")]] <- results$Coefficients
      all_plots[[paste(cond, rep, sep = "_")]] <- results$plot
      rep_ids <- c(rep_ids, paste(cond, rep, sep = "_"))
    }
  }

  message("\nFinalizing results...")
  
  # --------------------------------------------------------------------------
  # Reorder results to match original sample order
  # --------------------------------------------------------------------------
  NewData <- NewDataComplete[, Annotation$Sample, drop = FALSE]
  Propestimates <- PropestimatesComplete[, Annotation$Sample, drop = FALSE]
  
  # --------------------------------------------------------------------------
  # Compile fraction results into a single data frame
  # --------------------------------------------------------------------------
  Fraction <- do.call(rbind, FractionsComplete)
  Fraction <- as.data.frame(Fraction)
  # Calculate "Lost" fraction as remainder (1 - sum of observed fractions)
  Fraction$Lost <- rep(1, nrow(Fraction)) - rowSums(Fraction)
  Fraction$Replicate <- rownames(Fraction)
  rownames(Fraction) <- NULL
  
  # Compile coefficient results
  Coeff <- do.call(rbind, CoefficientComplete)
  Coeff <- as.data.frame(Coeff)
  Coeff$Replicate <- rownames(Coeff)
  rownames(Coeff) <- NULL

  message("FracFixR analysis complete!")
  
  # --------------------------------------------------------------------------
  # Return comprehensive results object
  # --------------------------------------------------------------------------
  list(
    OriginalData = DataNorm,
    Annotation = Annotation,
    Propestimates = Propestimates,
    NewData = NewData,
    Coefficients = Coeff,
    Fractions = Fraction,
    plots = all_plots
  )
}

# ==============================================================================
# SECTION 3: REPLICATE PROCESSING FUNCTION
# ==============================================================================

#' ProcessReplicate: Core NNLS Regression for Individual Replicates
#'
#' @description
#' This function implements the mathematical core of FracFixR: fitting a
#' non-negative least squares (NNLS) regression to estimate fraction weights
#' and correct individual transcript abundances.
#'
#' @param RepMat Data frame with transcripts as rows, samples as columns
#'   Must include a "Total" column representing the whole cell lysate
#' @param transcriptlist Character vector of transcript IDs to use for regression
#'   These should be informative transcripts in the appropriate abundance range
#'
#' @return List containing:
#'   - Propestimates: Proportion estimates for each transcript
#'   - NewData: Corrected count data
#'   - Coefficients: NNLS regression coefficients (fraction weights)
#'   - Fractions: Normalized fraction proportions
#'   - plot: Diagnostic plot of fitted vs residuals
#'
#' @details
#' Mathematical basis:
#' Total = α₀ + α₁×Fraction1 + α₂×Fraction2 + ... + ε
#' Where α₀ represents the "lost" fraction and other αᵢ are fraction weights
ProcessReplicate <- function(RepMat, transcriptlist) {
  Data <- data.frame(RepMat)
  n <- ncol(Data)

  # Validate presence of Total column (required for regression)
  if (!"Total" %in% colnames(RepMat)) {
    stop("RepMat must include a column named 'Total'. At least one sample of type 'Total' is required.")
  }

  # --------------------------------------------------------------------------
  # Prepare training data: subset to informative transcripts
  # --------------------------------------------------------------------------
  DataT <- Data[transcriptlist, ]
  # Remove any transcripts with NA values to ensure clean regression
  DataT <- DataT[rowSums(is.na(DataT)) == 0, ]
  
  # Create predictor matrix (all fractions except Total)
  X <- as.matrix(DataT[!is.element(colnames(DataT), "Total")])

  # --------------------------------------------------------------------------
  # Fit NNLS regression with intercept
  # The intercept captures the "lost" fraction not represented in sequenced samples
  # --------------------------------------------------------------------------
  X_with_intercept <- cbind(1, X)  # Add column of 1s for intercept
  fit <- nnls(X_with_intercept, DataT$Total)
  coef <- coef(fit)
  
  # Create diagnostic plot data
  PlotFit <- data.frame(
    fitted = as.vector(X_with_intercept %*% coef), 
    residuals = DataT$Total - as.vector(X_with_intercept %*% coef)
  )
  plotfit <- ggplot(PlotFit, aes(x = fitted, y = residuals)) + 
    geom_point() + 
    geom_hline(yintercept = 0)

  # --------------------------------------------------------------------------
  # Apply fitted model to all transcripts (not just training set)
  # --------------------------------------------------------------------------
  X_new <- as.matrix(Data[!is.element(colnames(Data), "Total")])
  seen_predict <- X_new %*% coef[-1]  # Predictions without intercept
  X_new <- cbind(1, X_new)
  predictions <- X_new %*% coef  # Full predictions including intercept

  # --------------------------------------------------------------------------
  # Calculate library sizes and fraction proportions
  # --------------------------------------------------------------------------
  Nlib <- colSums(Data)  # Library sizes for each sample
  Fractions <- c()
  Coefficients <- c(coef[1])  # Start with intercept (lost fraction)

  # Convert regression coefficients to fraction proportions
  # Account for different library sizes between Total and fractions
  for (j in 2:n) {
    # Fraction proportion = coefficient × fraction_library_size / total_library_size
    Fractions <- c(Fractions, coef[j] * Nlib[j] / Nlib[1])
    Coefficients <- c(Coefficients, coef[j])
  }
  names(Coefficients) <- colnames(Data)

  # --------------------------------------------------------------------------
  # Calculate proportion estimates for each transcript
  # --------------------------------------------------------------------------
  Propestimates <- Data
  Propestimates[, 1] <- ceiling(predictions)  # Total column gets predicted values
  
  # For each fraction, calculate proportion relative to predicted total
  # Use max() to ensure non-zero denominator
  for (j in 2:n) {
    Propestimates[, j] <- Data[, j] * coef[j] / apply(cbind(Data[, 1], seen_predict, 1), 1, max)
  }

  # --------------------------------------------------------------------------
  # Generate corrected count data
  # Normalize to match original Total library size
  # --------------------------------------------------------------------------
  NewData <- Propestimates
  for (j in 2:n) {
    # Scale fraction counts to match Total library size
    NewData[, j] <- ceiling(NewData[, j] * Nlib[1] / sum(NewData[, j]))
  }
  NewData[, 1] <- Data[, 1]  # Keep original Total counts

  # Return comprehensive results
  return(list(
    Propestimates = Propestimates,
    NewData = NewData,
    Coefficients = Coefficients,
    Fractions = Fractions,
    plot = plotfit
  ))
}

# ==============================================================================
# SECTION 4: DIFFERENTIAL PROPORTION TESTING
# ==============================================================================

#' DiffPropTest: Statistical Testing for Differential Proportions
#'
#' @description
#' Performs statistical testing to identify transcripts with significantly
#' different proportions between conditions in specified fraction(s).
#' Implements three test options: GLM (most powerful), Logit, and Wald.
#'
#' @param NormObject Output from FracFixR() function
#' @param Conditions Character vector of exactly 2 conditions to compare
#' @param Types Character vector of fraction type(s) to analyze
#'   Can be single fraction or multiple (will be combined)
#' @param Test Statistical test to use: "GLM", "Logit", or "Wald"
#'
#' @return Data frame with columns:
#'   - transcript: transcript identifier
#'   - mean_success_cond1/2: mean proportions in each condition
#'   - mean_diff: difference in proportions
#'   - log2FC: log2 fold change
#'   - pval: p-value from statistical test
#'   - padj: FDR-adjusted p-value
#'
#' @details
#' GLM: Uses binomial generalized linear model (most statistically powerful)
#' Logit: Faster alternative using logit transformation
#' Wald: Beta-binomial Wald test for overdispersed count data
DiffPropTest <- function(NormObject, Conditions, Types, Test = c("GLM", "Logit", "Wald")) {
  Test <- match.arg(Test)

  # Validate that Total samples exist in annotation
  if (!any(NormObject$Annotation$Type == "Total")) {
    stop("NormObject$Annotation must contain at least one sample of Type 'Total'.")
  }

  # Validate requested conditions exist
  if (!all(Conditions %in% unique(NormObject$Annotation$Condition))) {
    stop("Some specified Conditions are not found in NormObject$Annotation.")
  }

  # Validate requested types exist
  if (!all(Types %in% unique(NormObject$Annotation$Type))) {
    stop("Some specified Types are not found in NormObject$Annotation.")
  }
  
  # Create descriptive type label for messages
  if (length(Types)==1) type=Types else type=paste(Types[1], Types[2],sep="+")

  # Extract relevant data for the specified conditions and types
  Extraction <- extract_condition_matrix(NormObject$OriginalData, 
                                       NormObject$Propestimates, 
                                       NormObject$Annotation, 
                                       Conditions, 
                                       Types)

  # Prepare sample information with condition as factor
  sample_info <- Extraction$annotation %>%
    dplyr::select(Sample, Condition) %>%
    dplyr::mutate(Condition = factor(Condition, levels = Conditions))

  # --------------------------------------------------------------------------
  # Execute selected statistical test
  # --------------------------------------------------------------------------
  
  if (Test == "GLM") {
    message(paste('Performing GLM test comparing Conditions', Conditions[1], 'and', 
                Conditions[2], 'in the', type, 'fraction'))
    
    # Select transcripts with non-zero counts in all samples
    transcripts <- rownames(Extraction$counts[rowProds(as.matrix(Extraction$counts)) > 0, ])
    
    # Progress message for parallel processing
    message(sprintf("Testing %d transcripts using parallel processing...", length(transcripts)))
    
    # Run GLM in parallel for efficiency
    results_list <- future_lapply(transcripts, run_glm, 
                                counts = Extraction$counts, 
                                successes = Extraction$successes, 
                                sample_info = sample_info)
    results_df <- bind_rows(results_list)
  }

  if (Test == "Logit") {
    message(paste('Performing Logit test comparing Conditions', Conditions[1], 'and', 
                Conditions[2], 'in the', Types, 'fraction'))
    results_df <- logit_diff_test(counts = Extraction$counts, 
                                successes = Extraction$successes, 
                                annotation = sample_info, 
                                cond1 = Conditions[1], 
                                cond2 = Conditions[2])
  }

  if (Test == "Wald") {
    message(paste('Performing Wald test comparing Conditions', Conditions[1], 'and', 
                Conditions[2], 'in the', Types, 'fraction'))
    # Convert proportions to counts for beta-binomial model
    success_counts <- round(Extraction$successes * Extraction$counts)
    results_df <- beta_binomial_wald(counts = Extraction$counts, 
                                   successes = success_counts, 
                                   annotation = sample_info)
  }

  # Apply multiple testing correction (Benjamini-Hochberg FDR)
  message("Applying FDR correction...")
  results_df$padj <- p.adjust(results_df$pval, method = "BH")
  
  # Summary statistics
  n_sig <- sum(results_df$padj < 0.01, na.rm = TRUE)
  message(sprintf("Analysis complete. Found %d transcripts with padj < 0.01", n_sig))
  
  return(results_df)
}

# ==============================================================================
# SECTION 5: DATA EXTRACTION HELPER FUNCTION
# ==============================================================================

#' extract_condition_matrix: Extract and Prepare Data for Statistical Testing
#'
#' @description
#' Extracts count and proportion data for specified conditions and fraction types.
#' Handles both single fraction and combined fraction analyses.
#'
#' @param originalcounts Original count matrix
#' @param proportions Proportion estimates from FracFixR
#' @param annotation Sample annotation data frame
#' @param conditions Vector of conditions to extract
#' @param types Vector of fraction types to analyze
#'
#' @return List containing:
#'   - counts: Total counts from whole cell samples
#'   - successes: Proportion data for specified fractions
#'   - annotation: Filtered and processed annotation
extract_condition_matrix <- function(originalcounts, proportions, annotation, conditions, types) {
  # --------------------------------------------------------------------------
  # Extract Total samples for the specified conditions
  # These serve as the denominator in proportion calculations
  # --------------------------------------------------------------------------
  ann_wct <- annotation %>%
    dplyr::filter(Condition %in% conditions, Type == "Total")

  if (nrow(ann_wct) == 0) {
    stop("No samples with Type 'Total' found for the specified conditions.")
  }

  # Create count matrix from Total samples
  wct_samples <- ann_wct$Sample
  count_matrix <- originalcounts[, wct_samples, drop = FALSE]
  # Rename columns to include condition and replicate info
  colnames(count_matrix) <- paste(ann_wct$Condition, ann_wct$Replicate, sep = "_")

  # Filter annotation to specified conditions and types
  ann_filtered <- annotation %>%
    dplyr::filter(Condition %in% conditions, Type %in% types)

  # --------------------------------------------------------------------------
  # Handle single fraction vs combined fraction analysis
  # --------------------------------------------------------------------------
  
  if (length(types) == 1) {
    # Single fraction analysis: straightforward extraction
    selected_samples <- ann_filtered$Sample
    result_matrix <- proportions[, selected_samples, drop = FALSE]
    result_annotation <- ann_filtered %>%
      dplyr::filter(Sample %in% selected_samples)
    
    # Standardize column names
    colnames(result_matrix) <- paste(ann_filtered$Condition, ann_filtered$Replicate, sep = "_")
    result_annotation$Sample <- paste(ann_filtered$Condition, ann_filtered$Replicate, sep = "_")
    
  } else {
    # Multiple fraction analysis: sum proportions across fraction types
    result_list <- list()
    ann_list <- list()

    # Process each condition-replicate combination
    for (cond in unique(ann_filtered$Condition)) {
      for (rep in sort(unique(ann_filtered$Replicate))) {
        temp <- ann_filtered %>%
          dplyr::filter(Condition == cond, Replicate == rep)

        # Sum proportions across all requested fraction types
        summed_vector <- NULL
        for (t in types) {
          samples <- temp$Sample[temp$Type == t]
          if (length(samples) > 0) {
            vec <- rowSums(proportions[, samples, drop = FALSE])
            summed_vector <- if (is.null(summed_vector)) vec else summed_vector + vec
          }
        }

        # Store results with standardized naming
        colname <- paste(cond, rep, sep = "_")
        result_list[[colname]] <- summed_vector

        # Create corresponding annotation entry
        ann_list[[colname]] <- data.frame(
          Sample = colname,
          Condition = cond,
          Type = paste(types, collapse = "+"),
          Replicate = rep,
          stringsAsFactors = FALSE
        )
      }
    }

    # Combine results
    result_matrix <- do.call(cbind, result_list)
    result_annotation <- do.call(rbind, ann_list)
  }

  return(list(
    counts = count_matrix,
    successes = result_matrix,
    annotation = result_annotation
  ))
}

# ==============================================================================
# SECTION 6: STATISTICAL TEST IMPLEMENTATIONS
# ==============================================================================

#' run_glm: Binomial GLM for Single Transcript
#'
#' @description
#' Fits a binomial generalized linear model to test for differential
#' proportions between conditions for a single transcript.
#'
#' @param transcript Transcript identifier
#' @param counts Total count matrix
#' @param successes Proportion matrix
#' @param sample_info Sample metadata with Condition factor
#'
#' @return Data frame with test results for this transcript
run_glm <- function(transcript, counts, successes, sample_info) {
  # Validate transcript exists in data
  if (!(transcript %in% rownames(counts))) {
    stop(paste("Transcript", transcript, "not found in count matrix."))
  }

  # Ensure matrix alignment
  if (!all(colnames(counts) == colnames(successes))) {
    stop("Column names of 'counts' and 'successes' must match.")
  }

  # Verify all samples are in metadata
  if (!all(colnames(counts) %in% sample_info$Sample)) {
    stop("All sample names in counts/successes must be present in sample_info$Sample.")
  }

  # Extract data for this transcript
  valid_samples <- colnames(counts)
  total <- counts[transcript, valid_samples]
  prop_success <- successes[transcript, valid_samples]
  # Convert proportions to counts for binomial model
  success <- round(total * prop_success)
  failure <- total - success

  # Build data frame for GLM
  df_tmp <- data.frame(
    Sample = valid_samples,
    successes = as.numeric(success),
    failures = as.numeric(failure),
    prop = as.numeric(prop_success)
  ) %>%
    left_join(sample_info, by = "Sample")

  # Handle edge case: all zeros
  if (any(df_tmp$successes + df_tmp$failures == 0)) {
    return(data.frame(
      transcript = transcript,
      mean_success_cond1 = NA,
      mean_success_cond2 = NA,
      mean_diff = NA,
      log2FC = NA,
      pval = NA
    ))
  }

  # Fit GLM with error handling
  result <- tryCatch({
    # Binomial GLM: success/failure ~ Condition
    fit <- glm(cbind(successes, failures) ~ Condition, data = df_tmp, family = binomial)
    coeffs <- summary(fit)$coefficients
    # Extract p-value for condition effect
    pval <- coeffs[grep("^Condition", rownames(coeffs)), "Pr(>|z|)"][1]

    # Calculate group means
    group_means <- df_tmp %>%
      dplyr::group_by(Condition) %>%
      dplyr::summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop")

    # Ensure we have exactly 2 conditions
    if (nrow(group_means) != 2) {
      return(data.frame(
        transcript = transcript,
        mean_success_cond1 = NA,
        mean_success_cond2 = NA,
        mean_diff = NA,
        log2FC = NA,
        pval = NA
      ))
    }

    # Calculate effect sizes
    group_means <- group_means[order(group_means$Condition), ]
    mean_diff <- diff(group_means$mean_prop)
    log2FC <- log2(group_means$mean_prop[2] / group_means$mean_prop[1])

    data.frame(
      transcript = transcript,
      mean_success_cond1 = group_means$mean_prop[1],
      mean_success_cond2 = group_means$mean_prop[2],
      mean_diff = mean_diff,
      log2FC = log2FC,
      pval = pval
    )
  }, error = function(e) {
    # Return NA results if model fitting fails
    data.frame(
      transcript = transcript,
      mean_success_cond1 = NA,
      mean_success_cond2 = NA,
      mean_diff = NA,
      log2FC = NA,
      pval = NA
    )
  })

  return(result)
}

#' logit_diff_test: Logit-based Differential Test
#'
#' @description
#' Alternative to GLM using logit transformation. Faster but potentially
#' less powerful than the full GLM approach.
#'
#' @param counts Total count matrix
#' @param successes Proportion matrix  
#' @param annotation Sample metadata
#' @param cond1 First condition name
#' @param cond2 Second condition name
#'
#' @return Data frame with test results for all transcripts
logit_diff_test <- function(counts, successes, annotation, cond1, cond2) {
  # Validate matrix alignment
  if (!all(rownames(successes) == rownames(counts))) {
    stop("Row names of 'successes' and 'counts' must match.")
  }

  if (!all(colnames(successes) %in% annotation$Sample)) {
    stop("Not all column names in 'successes' are found in 'annotation$Sample'.")
  }

  # Process each transcript
  results <- lapply(rownames(counts), function(transcript) {
    # Build data frame for this transcript
    df <- data.frame(
      counts = as.numeric(counts[transcript, ]),
      successes = as.numeric(successes[transcript, ]),
      sample = colnames(counts)
    ) %>%
      dplyr::left_join(annotation, by = c("sample" = "Sample")) %>%
      dplyr::filter(Condition %in% c(cond1, cond2)) %>%
      dplyr::mutate(Condition = factor(Condition, levels = c(cond1, cond2)))

    # Check for valid comparison
    if (length(unique(df$Condition)) < 2) {
      return(data.frame(
        transcript = transcript,
        mean_success_cond1 = NA,
        mean_success_cond2 = NA,
        mean_diff = NA,
        log2FC = NA,
        pval = NA
      ))
    }

    # Fit model with error handling
    tryCatch({
      # Calculate condition means
      mean_success_cond1 <- mean(df$successes[df$Condition == cond1])
      mean_success_cond2 <- mean(df$successes[df$Condition == cond2])
      mean_diff <- mean_success_cond2 - mean_success_cond1

      # Fit binomial GLM
      model <- glm(cbind(successes, counts - successes) ~ Condition, 
                   family = binomial, data = df)
      summary_model <- summary(model)
      
      # Extract log2 fold change and p-value
      log2FC <- log2(exp(summary_model$coefficients[2, 1]))
      pval <- summary_model$coefficients[2, 4]

      data.frame(
        transcript = transcript,
        mean_success_cond1 = mean_success_cond1,
        mean_success_cond2 = mean_success_cond2,
        mean_diff = mean_diff,
        log2FC = log2FC,
        pval = pval
      )
    }, error = function(e) {
      data.frame(
        transcript = transcript,
        mean_success_cond1 = NA,
        mean_success_cond2 = NA,
        mean_diff = NA,
        log2FC = NA,
        pval = NA
      )
    })
  })

  do.call(rbind, results)
}

#' beta_binomial_wald: Beta-Binomial Wald Test
#'
#' @description
#' Implements Wald test using beta-binomial distribution to account for
#' overdispersion in count data. Useful when variance exceeds that expected
#' under binomial distribution.
#'
#' @param counts Total count matrix
#' @param successes Success count matrix (not proportions)
#' @param annotation Sample metadata
#'
#' @return Data frame with test results for all transcripts
beta_binomial_wald <- function(counts, successes, annotation) {
  # Validate matrix alignment
  if (!all(rownames(successes) == rownames(counts))) {
    stop("Row names of 'successes' and 'counts' must match.")
  }

  if (!all(colnames(successes) %in% annotation$Sample)) {
    stop("Not all samples in 'successes' matrix are found in 'annotation$Sample'.")
  }

  # Process each transcript
  results <- lapply(rownames(counts), function(transcript) {
    # Build data frame for this transcript
    df <- data.frame(
      counts = as.numeric(counts[transcript, ]),
      successes = as.numeric(successes[transcript, ]),
      sample = colnames(counts)
    ) %>%
      dplyr::left_join(annotation, by = c("sample" = "Sample"))

    # Ensure Condition is a factor
    df$Condition <- factor(df$Condition)
    cond_levels <- levels(df$Condition)

    # Check for valid comparison (need exactly 2 conditions)
    if (length(cond_levels) != 2) {
      return(data.frame(
        transcript = transcript,
        mean_success_cond1 = NA,
        mean_success_cond2 = NA,
        mean_diff = NA,
        log2FC = NA,
        pval = NA
      ))
    }

    # Fit beta-binomial model with error handling
    tryCatch({
      # Calculate condition means
      mean_success_cond1 <- mean(df$successes[df$Condition == cond_levels[1]])
      mean_success_cond2 <- mean(df$successes[df$Condition == cond_levels[2]])
      mean_diff <- mean_success_cond2 - mean_success_cond1

      # Fit beta-binomial model
      model <- aod::betabin(cbind(successes, counts - successes) ~ Condition, ~1, data = df)
      
      # Extract p-value and effect size (fixed string concatenation)
      coef_name <- paste0("Condition", cond_levels[2])
      pval <- summary(model)$coef[coef_name, "Pr(>Chi)"]
      log2FC <- log2(exp(coef(model)[coef_name]))

      data.frame(
        transcript = transcript,
        mean_success_cond1 = mean_success_cond1,
        mean_success_cond2 = mean_success_cond2,
        mean_diff = mean_diff,
        log2FC = log2FC,
        pval = pval
      )
    }, error = function(e) {
      data.frame(
        transcript = transcript,
        mean_success_cond1 = NA,
        mean_success_cond2 = NA,
        mean_diff = NA,
        log2FC = NA,
        pval = NA
      )
    })
  })

  do.call(rbind, results)
}

# ==============================================================================
# SECTION 7: VISUALIZATION FUNCTIONS
# ==============================================================================

#' PlotFractions: Visualize Fraction Proportions
#'
#' @description
#' Creates a stacked bar plot showing the distribution of RNA across fractions
#' for each replicate, including the "lost" fraction.
#'
#' @param FracFixed Output from FracFixR() function
#'
#' @return ggplot2 object showing fraction proportions
PlotFractions <- function(FracFixed) {
  df_temp <- FracFixed$Fractions
  
  # --------------------------------------------------------------------------
  # Reshape data from wide to long format for ggplot
  # --------------------------------------------------------------------------
  df_long <- df_temp %>%
    pivot_longer(
      cols = -Replicate,  # All columns except Replicate are proportions
      names_to = "Condition",
      values_to = "Proportion"
    )

  # --------------------------------------------------------------------------
  # Set up color scheme with "Lost" fraction always in grey
  # --------------------------------------------------------------------------
  conditions <- unique(df_long$Condition)
  other_conditions <- setdiff(conditions, "Lost")

  # Use colorblind-friendly palette for non-Lost conditions
  n_colors <- length(other_conditions)
  palette_colors <- brewer.pal(max(3, min(8, n_colors)), "Set2")[1:n_colors]
  names(palette_colors) <- other_conditions

  # Add grey for "Lost" fraction
  colors_named <- c(palette_colors, Lost = "grey80")

  # Set factor levels to ensure "Lost" appears at bottom of stack
  df_long <- df_long %>%
    mutate(Condition = factor(Condition, levels = c("Lost", other_conditions)))

  # --------------------------------------------------------------------------
  # Create stacked bar plot
  # --------------------------------------------------------------------------
  G <- ggplot(df_long, aes(x = Replicate, y = Proportion, fill = Condition)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(
      title = "Evaluation of Fraction Proportions",
      y = "Proportion",
      x = "Replicate"
    ) +
    scale_fill_manual(values = colors_named) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      text = element_text(size = 35)
    )

  return(G)
}

#' PlotComparison: Create Volcano Plot for Differential Results
#'
#' @description
#' Generates an enhanced volcano plot showing transcripts with significant
#' differential proportions between conditions.
#'
#' @param DiffPropResult Output from DiffPropTest() function
#' @param Conditions Character vector of conditions being compared
#' @param Types Character vector of fraction types analyzed
#' @param cutoff Optional y-axis maximum for plot
#'
#' @return EnhancedVolcano plot object
PlotComparison <- function(DiffPropResult, Conditions=NULL, Types=NULL, cutoff=NULL) {
  # --------------------------------------------------------------------------
  # Generate informative title and subtitle
  # --------------------------------------------------------------------------
  if (length(Types)==1) {
    subtitle <- paste(Types, 'comparison')
  } else if (length(Types)==2) {
    subtitle <- paste(paste(Types[1], Types[2], sep="+"), 'comparison')
  } else {
    subtitle <- NULL
  }
  
  if (length(Conditions)==2) {
    title <- paste(Conditions[1], 'vs', Conditions[2])
  } else {
    title <- paste('Comparison of fraction origin')
  }

  # --------------------------------------------------------------------------
  # Handle edge case of zero p-values
  # --------------------------------------------------------------------------
  if (min(DiffPropResult$padj, na.rm = TRUE) == 0) {
    warning(paste("One or more p-values is 0.", "Converting to 10^-1 * current", 
                  "lowest non-zero p-value..."), call. = FALSE)
    # Replace zeros with small non-zero value
    DiffPropResult$padj[which(DiffPropResult$padj == 0)] <- 
      min(DiffPropResult$padj[which(DiffPropResult$padj != 0)], na.rm = TRUE) * 10^-1
  }
  
  # Set y-axis limit
  if (length(cutoff) > 0) {
    ymax <- cutoff
  } else {
    ymax <- max(-log10(DiffPropResult$padj), na.rm = TRUE) + 5
  }
  
  # --------------------------------------------------------------------------
  # Create enhanced volcano plot
  # --------------------------------------------------------------------------
  PlotTrans <- EnhancedVolcano(DiffPropResult,
    lab = DiffPropResult$transcript,
    title = title,
    subtitle = subtitle,
    x = 'mean_diff',
    xlab = "Proportion shift",
    ylab = bquote(~-Log[10] ~ italic(P.adj)),
    legendLabels = c("NS", "prop shift", "adj.pval", "adj.pval + prop shift"),
    FCcutoff = 0.1,  # 10% proportion shift threshold
    xlim = c(min(DiffPropResult$mean_diff, na.rm = TRUE) - 0.1, 
             max(DiffPropResult$mean_diff, na.rm = TRUE) + 0.1),
    ylim = c(0, ymax),
    titleLabSize = 25,
    subtitleLabSize = 25,
    captionLabSize = 22,
    axisLabSize = 22,
    labSize = 0,  # Don't label individual points
    legendLabSize = 15,
    pCutoff = 0.01,  # Adjusted p-value threshold
    y = 'padj'  # Use adjusted p-values
  )
  
  return(PlotTrans)    
}