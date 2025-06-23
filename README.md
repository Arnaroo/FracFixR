# FracFixR: Fraction Correction Framework for RNA-seq Data

## Table of Contents

- [Introduction](#introduction)
- [The Problem](#the-problem)
- [How FracFixR Works](#how-fracfixr-works)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Quick Start](#quick-start)
- [Input Requirements](#input-requirements)
- [Output Description](#output-description)
- [Detailed Workflow](#detailed-workflow)
- [Function Reference](#function-reference)
- [Test Data Generation](#test-data-generation)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Contact](#contact)

## Introduction

FracFixR is a compositional statistical framework for absolute proportion estimation between fractions in RNA sequencing data. It addresses the fundamental challenge in fractionated RNA-seq experiments where library preparation and sequencing depth obscure the original proportions of RNA fractions.

### Key Features

- **Reconstructs original fraction proportions** using non-negative linear regression
- **Estimates the "lost" unrecoverable fraction** not captured in sequencing
- **Corrects individual transcript frequencies** for accurate abundance estimation
- **Performs differential proportion testing** between conditions
- **Handles any RNA fractionation protocol** (polysome profiling, subcellular localization, RNA-protein complexes)

### Applications

- Polysome profiling experiments
- Subcellular RNA localization studies
- Nuclear-cytoplasmic fractionation
- RNA-protein complex isolation
- Any experiment involving RNA fractionation

## The Problem

In fractionated RNA-seq experiments, several biases compromise data interpretation:

1. **Library preparation effects**: Amplification and sampling introduce variable efficiency across fractions
2. **Sequencing depth differences**: Different fractions may be sequenced at different depths
3. **Lost material**: Some RNA is unrecoverable during fractionation (e.g., degraded RNA, non-fractionated material)
4. **Compositional nature**: Fractions are parts of a whole, requiring specialized statistical treatment

Standard RNA-seq analysis tools (DESeq2, edgeR) are inappropriate because they:

- Assume independent samples rather than compositional data
- Cannot estimate the unobserved "lost" fraction
- Fail to account for global shifts in RNA distribution between conditions

## How FracFixR Works

FracFixR uses a mathematical model based on the compositional relationship between whole cell RNA and its fractions:

*Total = α₀ + α₁×Fraction1 + α₂×Fraction2 + ... + ε*

Where:
- *α₀* represents the "lost" fraction (intercept)
- *α₁, α₂, ...* are fraction weights
- *ε* is the error term

### The Algorithm

1. **Transcript Selection**: Identifies informative transcripts in the 70-96% abundance quantile range
2. **NNLS Regression**: Fits non-negative least squares model to estimate fraction weights
3. **Proportion Estimation**: Calculates individual transcript proportions in each fraction
4. **Count Correction**: Normalizes counts to account for library size differences
5. **Statistical Testing**: Performs differential proportion analysis using GLM/Logit/Wald tests

## Installation

### From GitHub (Recommended)

```r
install.packages("https://github.com/Arnaroo/FracFixR/releases/download/v1.0.0/FracFixR_1.0.0.tar.gz",repos = NULL,type = "source")```

### From GitHub (Development Version)

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install FracFixR from GitHub
devtools::install_github("Arnaroo/FracFixR")
```

### Verify Installation

```r
library(FracFixR)
packageVersion("FracFixR")
```

## Dependencies

FracFixR requires the following R packages:

### Core Dependencies
- **nnls** (≥ 1.4): Non-negative least squares regression
- **future.apply** (≥ 1.8.1): Parallel processing

### Data Manipulation
- **dplyr** (≥ 1.0.7): Data manipulation
- **tidyr** (≥ 1.1.3): Data tidying
- **matrixStats** (≥ 0.60.0): Efficient matrix operations

### Visualization
- **ggplot2** (≥ 3.3.5): Plotting framework
- **RColorBrewer** (≥ 1.1-2): Color palettes
- **EnhancedVolcano** (≥ 1.10.0): Volcano plots

### Statistical Testing
- **aod** (≥ 1.3.1): Beta-binomial regression

### System Requirements
- **R** ≥ 4.0.0
- **RAM**: Minimum 8GB (16GB recommended for large datasets)
- **CPU**: Multi-core processor recommended for parallel processing

## Quick Start

```r
library(FracFixR)

# Load your data
counts <- read.csv("counts_matrix.csv", row.names = 1)
annotation <- read.csv("annotation.csv")

# Convert to matrix
counts_matrix <- as.matrix(counts)

# Run FracFixR
results <- FracFixR(counts_matrix, annotation)

# Visualize fraction proportions
PlotFractions(results)

# Differential proportion testing
diff_results <- DiffPropTest(results, 
                            Conditions = c("Control", "Treatment"),
                            Types = "Heavy_Polysome",
                            Test = "GLM")

# Visualize differential results
PlotComparison(diff_results)
```

## Input Requirements

### 1. Count Matrix

A numeric matrix with:
- **Rows**: Transcripts/genes (must have unique row names)
- **Columns**: Samples (must match *Sample* column in annotation)
- **Values**: Raw read counts (not normalized)

```r
# Example structure
#         Sample1 Sample2 Sample3 Sample4
# Gene1      234     567     123     890
# Gene2      456     789     234     567
# Gene3      123     345     567     234
```

### 2. Annotation Data Frame

A data frame with required columns:
- **Sample**: Sample identifiers matching column names in count matrix
- **Condition**: Experimental condition (e.g., "Control", "Treatment")
- **Type**: Fraction type (must include at least one "Total")
- **Replicate**: Replicate identifier (e.g., "Rep1", "Rep2")

```r
# Example structure
#   Sample   Condition    Type     Replicate
# 1 Sample1  Control      Total    Rep1
# 2 Sample2  Control      Heavy_P  Rep1
# 3 Sample3  Treatment    Total    Rep1
# 4 Sample4  Treatment    Heavy_P  Rep1
```

### Important Notes
- Each condition-replicate combination must have one "Total" sample
- Fraction names must be consistent across replicates
- Minimum 2 replicates per condition recommended for statistical testing

## Output Description

### FracFixR Main Output

The `FracFixR()` function returns a list containing:

```r
results <- FracFixR(counts_matrix, annotation)

# Access components:
results$OriginalData    # Filtered input count matrix
results$Annotation      # Input annotation data
results$Propestimates   # Proportion estimates for each transcript
results$NewData         # Corrected count matrix
results$Coefficients    # Regression coefficients (α values)
results$Fractions       # Estimated fraction proportions
results$plots           # List of diagnostic plots
```

### DiffPropTest Output

The `DiffPropTest()` function returns a data frame with:

```r
# Column descriptions:
# transcript          - Transcript identifier
# mean_success_cond1  - Mean proportion in condition 1
# mean_success_cond2  - Mean proportion in condition 2
# mean_diff           - Difference in proportions (cond2 - cond1)
# log2FC              - Log2 fold change
# pval                - Raw p-value
# padj                - FDR-adjusted p-value
```

## Detailed Workflow

### Step 1: Prepare Your Data

```r
# Read count data
counts <- read.csv("polysome_counts.csv", row.names = 1)
counts_matrix <- as.matrix(counts)

# Create annotation
annotation <- data.frame(
  Sample = colnames(counts_matrix),
  Condition = c(rep("Control", 6), rep("Treatment", 6)),
  Type = rep(c("Total", "Light_Poly", "Heavy_Poly"), 4),
  Replicate = c(rep("Rep1", 3), rep("Rep2", 3), rep("Rep1", 3), rep("Rep2", 3))
)

# Verify data integrity
all(colnames(counts_matrix) %in% annotation$Sample)  # Should be TRUE
any(annotation$Type == "Total")                      # Should be TRUE
```

### Step 2: Run FracFixR

```r
# Run with default parameters
results <- FracFixR(MatrixCounts = counts_matrix, 
                   Annotation = annotation)

# The function will:
# 1. Filter transcripts present in Total samples
# 2. Process each condition-replicate combination
# 3. Fit NNLS regression
# 4. Calculate fraction proportions and lost fraction
# 5. Correct individual transcript abundances
```

### Step 3: Visualize Results

```r
# Plot fraction proportions
frac_plot <- PlotFractions(results)
ggsave("fraction_proportions.pdf", frac_plot, width = 10, height = 8)

# Access specific results
fraction_data <- results$Fractions
write.csv(fraction_data, "estimated_fractions.csv", row.names = FALSE)
```

### Step 4: Differential Testing

```r
# Compare conditions for specific fraction
diff_heavy <- DiffPropTest(NormObject = results,
                          Conditions = c("Control", "Treatment"),
                          Types = "Heavy_Poly",
                          Test = "GLM")

# Compare combined fractions
diff_combined <- DiffPropTest(NormObject = results,
                             Conditions = c("Control", "Treatment"),
                             Types = c("Light_Poly", "Heavy_Poly"),
                             Test = "GLM")

# Filter significant results
sig_transcripts <- diff_heavy[diff_heavy$padj < 0.01, ]
```

### Step 5: Generate Reports

```r
# Volcano plot
volcano <- PlotComparison(diff_heavy, 
                         Conditions = c("Control", "Treatment"),
                         Types = "Heavy_Poly",
                         cutoff = 50)

ggsave("volcano_plot.pdf", volcano, width = 12, height = 10)

# Export results
write.csv(diff_heavy, "differential_proportion_results.csv", row.names = FALSE)
```

## Function Reference

### Main Functions

#### 1. FracFixR()

```r
FracFixR(MatrixCounts, Annotation)
```

**Parameters:**
- `MatrixCounts`: Numeric matrix of raw counts (genes × samples)
- `Annotation`: Data frame with Sample, Condition, Type, and Replicate columns

**Returns:** List containing corrected data and analysis results

**Details:**
- Automatically detects and uses available CPU cores for parallel processing
- Filters transcripts to those present in Total samples
- Uses 70-96% quantile range for regression transcript selection

#### 2. DiffPropTest()

```r
DiffPropTest(NormObject, Conditions, Types, Test = c("GLM", "Logit", "Wald"))
```

**Parameters:**
- `NormObject`: Output from FracFixR()
- `Conditions`: Character vector of exactly 2 conditions to compare
- `Types`: Character vector of fraction type(s) to analyze
- `Test`: Statistical test method
  - `"GLM"`: Binomial generalized linear model (most powerful, default)
  - `"Logit"`: Logit transformation test (faster)
  - `"Wald"`: Beta-binomial Wald test (for overdispersed data)

**Returns:** Data frame with differential proportion results

### Visualization Functions

#### 3. PlotFractions()

```r
PlotFractions(FracFixed)
```

**Parameters:**
- `FracFixed`: Output from FracFixR()

**Returns:** ggplot2 object showing stacked bar plot of fraction proportions

#### 4. PlotComparison()

```r
PlotComparison(DiffPropResult, Conditions = NULL, Types = NULL, cutoff = NULL)
```

**Parameters:**
- `DiffPropResult`: Output from DiffPropTest()
- `Conditions`: Optional character vector for plot title
- `Types`: Optional character vector for plot subtitle
- `cutoff`: Optional y-axis maximum for volcano plot

**Returns:** Enhanced volcano plot

## Test Data Generation

FracFixR includes a Python script for generating synthetic test data with known ground truth:

### Prerequisites

```bash
pip install numpy pandas
```

### Generate Test Data

```bash
# generate_data.py usage
python generate_data.py \
  --total_reads 1000000 \
  --n_transcripts 10000 \
  --fraction_weights "0.3,0.5" \
  --lost_fraction 0.2 \
  --output_dir test_data/
```

### Parameters:
- `--total_reads`: Total number of reads in dataset
- `--n_transcripts`: Number of transcripts to simulate
- `--fraction_weights`: Comma-separated fraction proportions
- `--lost_fraction`: Proportion of lost/unrecoverable material
- `--recovery_rates`: Library recovery rates for each fraction
- `--output_dir`: Output directory for generated files

### Load Test Data in R

```r
# Load generated test data
test_counts <- read.csv("test_data/counts_matrix.csv", row.names = 1)
test_annotation <- read.csv("test_data/annotation.csv")
ground_truth <- read.csv("test_data/ground_truth.csv")

# Run FracFixR
test_results <- FracFixR(as.matrix(test_counts), test_annotation)

# Compare with ground truth
cor(test_results$Fractions$Lost[1], ground_truth$lost_fraction[1])
```

## Examples

### Example 1: Basic Polysome Profiling

```r
# Minimal example with polysome data
library(FracFixR)

# Create example data
set.seed(123)
n_genes <- 1000
n_samples <- 12

# Generate count matrix
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
  Condition = rep(c("Control", "Stressed"), each = 6),
  Type = rep(c("Total", "Monosome", "Polysome"), 4),
  Replicate = c(rep("Rep1", 3), rep("Rep2", 3), rep("Rep1", 3), rep("Rep2", 3))
)

# Run analysis
results <- FracFixR(counts, annotation)

# Check results
print(results$Fractions)
PlotFractions(results)
```

### Example 2: Subcellular Fractionation

```r
# Nuclear-cytoplasmic fractionation example
annotation_subcell <- data.frame(
  Sample = colnames(counts_matrix),
  Condition = rep(c("WT", "Mutant"), each = 6),
  Type = rep(c("Total", "Nuclear", "Cytoplasmic"), 4),
  Replicate = rep(c("Rep1", "Rep2"), each = 3, times = 2)
)

# Run FracFixR
subcell_results <- FracFixR(counts_matrix, annotation_subcell)

# Test for differential nuclear localization
diff_nuclear <- DiffPropTest(subcell_results,
                            Conditions = c("WT", "Mutant"),
                            Types = "Nuclear",
                            Test = "GLM")

# Find transcripts with altered localization
nuclear_shifted <- diff_nuclear[abs(diff_nuclear$mean_diff) > 0.2 & 
                               diff_nuclear$padj < 0.01, ]
```

### Example 3: Multiple Fraction Analysis

```r
# Combine multiple fractions for analysis
combined_poly <- DiffPropTest(results,
                             Conditions = c("Control", "Stressed"),
                             Types = c("Monosome", "Polysome"),
                             Test = "GLM")

# This tests for overall ribosome association changes
# regardless of whether transcripts shift between mono- and polysomes
```

## Troubleshooting

### Common Issues and Solutions

#### 1. "No transcripts selected for regression"

**Cause**: Too few transcripts in the 70-96% quantile range

**Solution**:
```r
# Check transcript abundance distribution
total_counts <- rowSums(counts_matrix[, annotation$Type == "Total"])
hist(log10(total_counts + 1), breaks = 50)

# Adjust if needed (not recommended unless necessary)
# Contact package maintainers for guidance
```

#### 2. "Column names must match annotation Sample"

**Cause**: Mismatch between count matrix columns and annotation

**Solution**:
```r
# Check for mismatches
setdiff(colnames(counts_matrix), annotation$Sample)
setdiff(annotation$Sample, colnames(counts_matrix))

# Ensure exact matching
colnames(counts_matrix) <- gsub(" ", "_", colnames(counts_matrix))
annotation$Sample <- gsub(" ", "_", annotation$Sample)
```

#### 3. Memory issues with large datasets

**Solution**:
```r
# Reduce parallel workers
library(future)
plan("multisession", workers = 2)  # Use only 2 cores

# Process in batches
# Contact maintainers for batch processing scripts
```

#### 4. Convergence warnings in GLM

**Cause**: Extreme proportions or low counts

**Solution**:
```r
# Try alternative test
diff_results <- DiffPropTest(results,
                            Conditions = c("Control", "Treatment"),
                            Types = "Heavy_Poly",
                            Test = "Wald")  # More robust to extremes
```

### Performance Tips

1. **Parallel Processing**: FracFixR automatically uses available cores
   ```r
   # Check current plan
   future::plan()
   
   # Adjust if needed
   future::plan("multisession", workers = 4)
   ```

2. **Memory Management**: For datasets >50,000 transcripts
   ```r
   # Pre-filter low-abundance transcripts
   keep <- rowSums(counts_matrix) >= 10
   counts_filtered <- counts_matrix[keep, ]
   ```

3. **Speed Optimization**: Use Logit test for initial exploration
   ```r
   # Fast initial scan
   quick_test <- DiffPropTest(results, 
                             Conditions = c("A", "B"),
                             Types = "Fraction1",
                             Test = "Logit")
   ```

## Citation

If you use FracFixR in your research, please cite:

> Cleynen A, Ravindran A, Shirokikh N (2024). FracFixR: A compositional 
> statistical framework for absolute proportion estimation between fractions 
> in RNA sequencing data. Bioinformatics, btx###. 
> doi: 10.1093/bioinformatics/btx###

## Contact

- **Bug Reports**: https://github.com/Arnaroo/FracFixR/issues
- **Feature Requests**: https://github.com/Arnaroo/FracFixR/issues
- **General Questions**: fracfixr@gmail.com

### Authors

- Alice Cleynen (alice.cleynen@umontpellier.fr) - Statistical methodology
- Agin Ravindran - Experimental data generation
- Nikolay Shirokikh (nikolay.shirokikh@uwa.edu.au) - Conceptual framework

### Contributing

We welcome contributions! Please see our Contributing Guidelines for details.

### License

FracFixR is released under the CC BY-NC ND License. See LICENSE file for details.
