################################################################################
# FracFixR Complete Workflow Example
################################################################################
# This script demonstrates a complete analysis workflow using FracFixR
# for polysome profiling data analysis
#
# Author: FracFixR Team
# Date: 2024
################################################################################

# Load required libraries
library(FracFixR)
library(ggplot2)
library(dplyr)

# Set working directory (adjust as needed)
# setwd("~/polysome_analysis")

# ==============================================================================
# STEP 1: Load and Prepare Data
# ==============================================================================

# Option A: Use example data included with package
data(example_counts)
data(example_annotation)

# Option B: Load your own data
# counts <- read.csv("counts_matrix.csv", row.names = 1)
# annotation <- read.csv("sample_annotation.csv")
# 
# # Convert to matrix if needed
# counts_matrix <- as.matrix(counts)

# Examine the data structure
cat("Data dimensions:\n")
cat(sprintf("  Genes: %d\n", nrow(example_counts)))
cat(sprintf("  Samples: %d\n", ncol(example_counts)))
cat("\nAnnotation summary:\n")
print(table(example_annotation$Condition, example_annotation$Type))

# ==============================================================================
# STEP 2: Quality Control
# ==============================================================================

# Check library sizes
lib_sizes <- colSums(example_counts)
cat("\nLibrary sizes:\n")
print(summary(lib_sizes))

# Visualize library sizes
pdf("library_sizes.pdf", width = 10, height = 6)
par(mar = c(8, 4, 2, 2))
barplot(lib_sizes, 
        las = 2, 
        main = "Library Sizes by Sample",
        ylab = "Total Counts",
        col = rainbow(length(lib_sizes)))
dev.off()

# Check for genes with zero counts in Total samples
total_samples <- example_annotation$Type == "Total"
zero_genes <- sum(rowSums(example_counts[, total_samples]) == 0)
cat(sprintf("\nGenes with zero counts in Total samples: %d\n", zero_genes))

# Filter low-count genes (optional)
min_count <- 10
keep <- rowSums(example_counts) >= min_count
filtered_counts <- example_counts[keep, ]
cat(sprintf("Genes after filtering: %d\n", nrow(filtered_counts)))

# ==============================================================================
# STEP 3: Run FracFixR
# ==============================================================================

cat("\nRunning FracFixR analysis...\n")
fracfixr_results <- FracFixR(
  MatrixCounts = example_counts,
  Annotation = example_annotation
)

# Examine the results structure
cat("\nFracFixR output components:\n")
print(names(fracfixr_results))

# ==============================================================================
# STEP 4: Visualize Fraction Proportions
# ==============================================================================

# Plot fraction proportions
frac_plot <- PlotFractions(fracfixr_results)

# Customize the plot
frac_plot_custom <- frac_plot +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold")
  ) +
  scale_fill_brewer(palette = "Set3")

# Save the plot
ggsave("fraction_proportions.pdf", frac_plot_custom, width = 10, height = 8)

# Extract fraction data for further analysis
fraction_data <- fracfixr_results$Fractions
cat("\nEstimated fraction proportions:\n")
print(fraction_data)

# Calculate mean proportions by condition
mean_fractions <- fraction_data %>%
  tidyr::separate(Replicate, into = c("Condition", "Rep"), sep = "_") %>%
  group_by(Condition) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

cat("\nMean fraction proportions by condition:\n")
print(mean_fractions)

# ==============================================================================
# STEP 5: Differential Proportion Testing
# ==============================================================================

# Test for differential heavy polysome association
cat("\nTesting for differential heavy polysome association...\n")
diff_heavy <- DiffPropTest(
  NormObject = fracfixr_results,
  Conditions = c("Control", "Treatment"),
  Types = "Heavy_Polysome",
  Test = "GLM"
)

# Summary of results
n_tested <- nrow(diff_heavy)
n_sig_01 <- sum(diff_heavy$padj < 0.01, na.rm = TRUE)
n_sig_05 <- sum(diff_heavy$padj < 0.05, na.rm = TRUE)

cat(sprintf("\nDifferential testing summary:\n"))
cat(sprintf("  Transcripts tested: %d\n", n_tested))
cat(sprintf("  Significant at FDR < 0.01: %d\n", n_sig_01))
cat(sprintf("  Significant at FDR < 0.05: %d\n", n_sig_05))

# Get top differentially associated transcripts
top_genes <- diff_heavy %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(20)

cat("\nTop 20 differentially associated transcripts:\n")
print(top_genes[, c("transcript", "mean_diff", "log2FC", "padj")])

# Save full results
write.csv(diff_heavy, "differential_heavy_polysome.csv", row.names = FALSE)

# ==============================================================================
# STEP 6: Create Volcano Plot
# ==============================================================================

# Basic volcano plot
volcano <- PlotComparison(
  diff_heavy,
  Conditions = c("Control", "Treatment"),
  Types = "Heavy_Polysome"
)

# Save volcano plot
ggsave("volcano_heavy_polysome.pdf", volcano, width = 12, height = 10)

# ==============================================================================
# STEP 7: Test Combined Fractions
# ==============================================================================

# Test overall ribosome association (light + heavy)
cat("\nTesting for differential overall ribosome association...\n")
diff_combined <- DiffPropTest(
  NormObject = fracfixr_results,
  Conditions = c("Control", "Treatment"),
  Types = c("Light_Polysome", "Heavy_Polysome"),
  Test = "GLM"
)

# Summary
n_sig_combined <- sum(diff_combined$padj < 0.01, na.rm = TRUE)
cat(sprintf("Transcripts with altered ribosome association: %d\n", n_sig_combined))

# Compare single vs combined fraction results
merged_results <- merge(
  diff_heavy[, c("transcript", "mean_diff", "padj")],
  diff_combined[, c("transcript", "mean_diff", "padj")],
  by = "transcript",
  suffixes = c("_heavy", "_combined")
)

# Identify transcripts significant in both analyses
both_sig <- merged_results %>%
  filter(padj_heavy < 0.01 & padj_combined < 0.01)

cat(sprintf("\nTranscripts significant in both analyses: %d\n", nrow(both_sig)))

# ==============================================================================
# STEP 8: Export Corrected Counts
# ==============================================================================

# Get corrected counts
corrected_counts <- fracfixr_results$NewData

# Save corrected counts
write.csv(corrected_counts, "corrected_counts.csv")

# Compare original vs corrected for a specific gene
gene_example <- "Gene1"
original <- example_counts[gene_example, ]
corrected <- corrected_counts[gene_example, ]

comparison <- data.frame(
  Sample = names(original),
  Original = as.numeric(original),
  Corrected = as.numeric(corrected),
  Type = example_annotation$Type
)

cat(sprintf("\nExample comparison for %s:\n", gene_example))
print(comparison)

# ==============================================================================
# STEP 9: Generate Summary Report
# ==============================================================================

# Create a summary report
report <- list(
  date = Sys.Date(),
  n_genes = nrow(example_counts),
  n_samples = ncol(example_counts),
  conditions = unique(example_annotation$Condition),
  fraction_types = unique(example_annotation$Type),
  mean_lost_fraction = mean(fraction_data$Lost),
  n_diff_heavy = n_sig_01,
  n_diff_combined = n_sig_combined
)

# Save summary
saveRDS(report, "analysis_summary.rds")

cat("\n===============================================\n")
cat("Analysis complete! Generated files:\n")
cat("  - library_sizes.pdf\n")
cat("  - fraction_proportions.pdf\n")
cat("  - differential_heavy_polysome.csv\n")
cat("  - volcano_heavy_polysome.pdf\n")
cat("  - corrected_counts.csv\n")
cat("  - analysis_summary.rds\n")
cat("===============================================\n")

# ==============================================================================
# STEP 10: Session Information
# ==============================================================================

cat("\nSession information:\n")
sessionInfo()