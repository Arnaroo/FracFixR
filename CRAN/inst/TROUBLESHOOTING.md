# FracFixR Troubleshooting Guide

## Table of Contents
1. [Installation Issues](#installation-issues)
2. [Data Input Errors](#data-input-errors)
3. [Runtime Errors](#runtime-errors)
4. [Statistical Testing Issues](#statistical-testing-issues)
5. [Memory and Performance](#memory-and-performance)
6. [Interpretation Questions](#interpretation-questions)

## Installation Issues

### Problem: Package won't install from CRAN
```
Error: package 'FracFixR' is not available
```

**Solution:**
```r
# Check your R version
R.version.string  # Should be >= 4.0.0

# Try different CRAN mirror
chooseCRANmirror()
install.packages("FracFixR")

# Or install from GitHub
devtools::install_github("Arnaroo/FracFixR")
```

### Problem: Dependencies fail to install
```
Error: dependency 'EnhancedVolcano' is not available
```

**Solution:**
```r
# EnhancedVolcano is from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")

# Then retry FracFixR installation
install.packages("FracFixR")
```

## Data Input Errors

### Problem: "All MatrixCounts columns must be in Annotation$Sample"

**Diagnosis:**
```r
# Check for mismatches
setdiff(colnames(counts), annotation$Sample)  # Shows columns not in annotation
setdiff(annotation$Sample, colnames(counts))  # Shows annotations without data
```

**Solutions:**
```r
# Fix whitespace issues
colnames(counts) <- trimws(colnames(counts))
annotation$Sample <- trimws(annotation$Sample)

# Fix case sensitivity
colnames(counts) <- gsub("sample", "Sample", colnames(counts))

# Check exact matches
all(colnames(counts) %in% annotation$Sample)  # Should be TRUE
```

### Problem: "Annotation must include at least one sample with Type == 'Total'"

**Solution:**
```r
# Check Type values
unique(annotation$Type)

# Common fixes:
# 1. Case sensitivity
annotation$Type <- gsub("total", "Total", annotation$Type)

# 2. Whitespace
annotation$Type <- trimws(annotation$Type)

# 3. Typos
annotation$Type[annotation$Type == "Totla"] <- "Total"

# Verify
any(annotation$Type == "Total")  # Should be TRUE
```

### Problem: "MatrixCounts must be a numeric matrix"

**Solution:**
```r
# Convert data.frame to matrix
counts <- as.matrix(counts)

# Ensure numeric
counts <- apply(counts, 2, as.numeric)

# Check for non-numeric values
which(is.na(counts), arr.ind = TRUE)  # Shows problematic cells

# Remove non-numeric rows/columns if needed
numeric_cols <- sapply(counts, is.numeric)
counts <- counts[, numeric_cols]
```

## Runtime Errors

### Problem: "No transcripts selected for regression in condition X"

**Diagnosis:**
```r
# Check transcript abundance distribution
total_samples <- annotation$Type == "Total" & annotation$Condition == "X"
total_counts <- rowSums(counts[, total_samples, drop = FALSE])
hist(log10(total_counts + 1), breaks = 50)
quantile(total_counts, c(0.7, 0.96))  # Default selection range
```

**Solutions:**
```r
# 1. Check if you have enough expressed genes
sum(total_counts > 0)  # Should be > 10

# 2. Check for sample swaps
# Ensure Total samples have higher counts than fractions

# 3. If legitimate low expression, contact maintainers for guidance
```

### Problem: "Each replicate must have at least two types"

**Solution:**
```r
# Check your annotation structure
table(annotation$Replicate, annotation$Type)

# Each replicate should have Total + at least one fraction
# Fix by ensuring complete experimental design
```

### Problem: Convergence warnings in GLM

**Solutions:**
```r
# 1. Try Wald test (more robust to extreme values)
diff_results <- DiffPropTest(results,
                            Conditions = c("A", "B"),
                            Types = "Fraction",
                            Test = "Wald")

# 2. Filter genes with very low counts
min_count_threshold <- 10
keep <- rowSums(counts) >= min_count_threshold
filtered_counts <- counts[keep, ]

# 3. Check for outlier samples
boxplot(log10(colSums(counts) + 1), las = 2)
```

## Statistical Testing Issues

### Problem: All p-values are NA

**Diagnosis:**
```r
# Check if you have replicates
table(annotation$Condition, annotation$Replicate)

# Check for zero counts in tested fractions
fraction_samples <- annotation$Type == "YourFraction"
sum(rowSums(counts[, fraction_samples]) > 0)
```

**Solution:**
- Ensure at least 2 replicates per condition
- Check that the tested fraction has non-zero counts
- Verify fraction names match exactly

### Problem: No significant results

**Checks:**
```r
# 1. Verify biological variability exists
# Plot PCA of your fractions
pca_data <- prcomp(t(log10(counts + 1)))
plot(pca_data$x[,1], pca_data$x[,2])

# 2. Check proportion distributions
hist(diff_results$mean_diff, breaks = 50)

# 3. Try less stringent threshold
sum(diff_results$pval < 0.05, na.rm = TRUE)  # Before multiple testing
```

### Problem: Different results with different tests

This is expected! The tests have different assumptions:
- **GLM**: Assumes binomial distribution, most powerful
- **Logit**: Uses transformation, faster but less powerful
- **Wald**: Accounts for overdispersion, more conservative

## Memory and Performance

### Problem: Out of memory errors

**Solutions:**
```r
# 1. Reduce parallel workers
library(future)
plan(multisession, workers = 2)  # Use only 2 cores

# 2. Pre-filter genes
keep <- rowSums(counts) >= 10  # Minimum total count
filtered_counts <- counts[keep, ]

# 3. Process in batches (contact maintainers for scripts)
```

### Problem: Analysis takes too long

**Solutions:**
```r
# 1. Check parallel processing is working
future::availableCores()

# 2. Use Logit test for initial exploration
diff_quick <- DiffPropTest(results, 
                          Conditions = c("A", "B"),
                          Types = "Fraction",
                          Test = "Logit")  # Faster

# 3. Reduce gene number for testing
test_counts <- counts[1:1000, ]  # Subset for testing
```

## Interpretation Questions

### Q: What does the "Lost" fraction represent?

The lost fraction includes:
- RNA lost during fractionation procedure
- RNA in fractions not sequenced (e.g., free RNA in polysome profiling)
- Technical losses during library preparation
- Degraded or damaged RNA

High lost fraction (>50%) may indicate technical issues.

### Q: How do I interpret mean_diff?

- **Positive mean_diff**: Higher proportion in condition 2
- **Negative mean_diff**: Higher proportion in condition 1
- **Magnitude**: Absolute change in proportion (0.1 = 10% change)

Example:
```
mean_diff = 0.15, padj < 0.01
→ Gene has 15% higher proportion in condition 2 (significant)
```

### Q: Can I use normalized counts as input?

**No!** FracFixR requires raw counts because:
- It performs its own compositional normalization
- Standard normalizations assume independence between samples
- Fractions are compositionally related (parts of a whole)

### Q: How many replicates do I need?

- **Minimum**: 2 per condition (will get warning)
- **Recommended**: 3+ per condition
- **For subtle changes**: 4-6 per condition

Power depends on:
- Effect size (proportion changes)
- Technical variability
- Number of genes tested

### Q: Can I compare more than 2 conditions?

Currently, DiffPropTest compares exactly 2 conditions at a time.

For multiple conditions:
```r
# Pairwise comparisons
diff_AvsB <- DiffPropTest(results, Conditions = c("A", "B"), ...)
diff_AvsC <- DiffPropTest(results, Conditions = c("A", "C"), ...)
diff_BvsC <- DiffPropTest(results, Conditions = c("B", "C"), ...)
```

## Getting Help

1. **Check documentation**:
   ```r
   ?FracFixR
   vignette("FracFixR-intro")
   ```

2. **Reproducible example**:
   ```r
   # Create minimal example
   set.seed(123)
   mini_counts <- matrix(rpois(100, 50), nrow = 10)
   # ... show the error
   ```

3. **Report issues**:
   - GitHub: https://github.com/Arnaroo/FracFixR/issues
   - Include: sessionInfo(), error message, minimal example

4. **Contact maintainers**:
   - General questions: fracfixr@gmail.com
   - Bug reports: Use GitHub issues

## Common Warnings (Usually Safe)

1. **"Some conditions have fewer than 2 replicates"**
   - Analysis will run but with limited statistical power

2. **"Only X transcripts selected for condition Y"**
   - If X > 10, analysis should be reliable
   - Check your data if X < 10

3. **"Converting 0 p-values"**
   - Very significant results, p-values below machine precision
   - Converted to small non-zero values for plotting