# FracFixR Quick Reference Guide

## Installation
```r
# From CRAN
install.packages("FracFixR")

# From GitHub
devtools::install_github("Arnaroo/FracFixR")
```

## Basic Workflow

### 1. Load Package and Data
```r
library(FracFixR)

# Your data
counts <- as.matrix(read.csv("counts.csv", row.names = 1))
annotation <- read.csv("annotation.csv")
```

### 2. Run FracFixR
```r
results <- FracFixR(MatrixCounts = counts, 
                    Annotation = annotation)
```

### 3. Visualize Fractions
```r
PlotFractions(results)
```

### 4. Differential Testing
```r
diff_results <- DiffPropTest(results,
                            Conditions = c("Control", "Treatment"),
                            Types = "Heavy_Polysome",
                            Test = "GLM")
```

### 5. Volcano Plot
```r
PlotComparison(diff_results,
               Conditions = c("Control", "Treatment"),
               Types = "Heavy_Polysome")
```

## Key Functions

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `FracFixR()` | Main analysis | `MatrixCounts`, `Annotation` |
| `DiffPropTest()` | Differential testing | `NormObject`, `Conditions`, `Types`, `Test` |
| `PlotFractions()` | Visualize proportions | `FracFixed` |
| `PlotComparison()` | Volcano plot | `DiffPropResult`, `Conditions`, `Types` |

## Annotation Requirements

Required columns:
- **Sample**: Must match column names in count matrix
- **Condition**: Experimental conditions
- **Type**: Fraction types (must include "Total")
- **Replicate**: Replicate identifiers

Example:
```
Sample    Condition  Type    Replicate
Sample1   Control    Total   Rep1
Sample2   Control    Light   Rep1
Sample3   Control    Heavy   Rep1
```

## Statistical Tests

- **GLM** (default): Most powerful, binomial GLM
- **Logit**: Faster, logit transformation
- **Wald**: Beta-binomial for overdispersed data

## Common Analyses

### Single Fraction
```r
diff_heavy <- DiffPropTest(results,
                          Conditions = c("A", "B"),
                          Types = "Heavy")
```

### Combined Fractions
```r
diff_combined <- DiffPropTest(results,
                             Conditions = c("A", "B"),
                             Types = c("Light", "Heavy"))
```

### Filter Results
```r
# Significant genes
sig_genes <- diff_results[diff_results$padj < 0.01, ]

# Top changed genes
top_genes <- diff_results[order(abs(diff_results$mean_diff), 
                                decreasing = TRUE), ][1:50, ]
```

## Output Components

`FracFixR()` returns:
- `$OriginalData`: Filtered input counts
- `$Annotation`: Sample annotation
- `$Propestimates`: Proportion estimates
- `$NewData`: Corrected counts
- `$Coefficients`: Regression coefficients
- `$Fractions`: Fraction proportions
- `$plots`: Diagnostic plots

`DiffPropTest()` returns:
- `transcript`: Gene/transcript ID
- `mean_success_cond1/2`: Mean proportions
- `mean_diff`: Difference in proportions
- `log2FC`: Log2 fold change
- `pval`: Raw p-value
- `padj`: FDR-adjusted p-value

## Tips

1. **Check Total samples**: Ensure each condition-replicate has a "Total"
2. **Minimum replicates**: Use ≥2 replicates per condition
3. **Parallel processing**: Automatic, uses available cores
4. **Memory**: For >50k genes, pre-filter low counts
5. **Interpretation**: Positive mean_diff = higher in condition 2

## Example Datasets

```r
# Polysome profiling
data(example_counts)
data(example_annotation)

# Alternative annotations
data(polysome_annotation)    # Monosome/Polysome
data(subcellular_annotation)  # Nuclear/Cytoplasmic
```

## Troubleshooting

See full troubleshooting guide in package documentation or:
```r
vignette("FracFixR-intro")
help(FracFixR)
```

## Citation

```r
citation("FracFixR")
```

Cleynen et al. (2024) FracFixR: A compositional statistical framework for absolute proportion estimation between fractions in RNA sequencing data. *Bioinformatics*.