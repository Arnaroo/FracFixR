# FracFixR 1.1.0

## Improvements and New Features

### Bug fixes
* Fixed parallel worker setup to use `parallelly::availableCores()` instead of
  `parallel::detectCores() - 1`. This avoids exceeding the hard localhost worker
  limit enforced by `parallelly` and respects the `mc.cores` option set by
  CRAN/CI environments.

### New features
* `FracFixR()` gains two new parameters `st1` (default 0.6) and `st2` (default
  0.999) that allow the user to control the quantile range used to select
  informative transcripts for the NNLS regression fit. Previously these thresholds
  were fixed at 70% and 96%.
* New exported function `get_corrected_counts()`: converts the proportion matrix
  from `FracFixR()` back to an interpretable count matrix by multiplying each
  sample's proportions by the Total abundance of the matched replicate.

### Internal changes
* `TotalSum` computation now uses `na.rm = TRUE` in `rowSums()` and transcript
  selection now explicitly filters out `NA` values, improving robustness when
  Total samples contain missing data.

---

# FracFixR 1.0.0

## Initial CRAN Release

This is the first release of FracFixR, a compositional statistical framework for absolute proportion estimation between fractions in RNA sequencing data.

### Features

* **Core functionality**
  - `FracFixR()`: Main function for fraction correction using NNLS regression
  - `DiffPropTest()`: Statistical testing for differential proportions between conditions
  - Support for any RNA fractionation protocol (polysome profiling, subcellular localization, etc.)

* **Statistical methods**
  - Non-negative least squares (NNLS) regression for fraction weight estimation
  - Three test options: GLM (binomial), Logit, and Beta-binomial Wald
  - Automatic estimation of "lost" unrecoverable fraction
  - FDR correction for multiple testing

* **Visualization**
  - `PlotFractions()`: Stacked bar plots of fraction proportions
  - `PlotComparison()`: Enhanced volcano plots for differential results
  - Diagnostic plots for regression quality

* **Performance**
  - Automatic parallel processing using available CPU cores
  - Efficient handling of large datasets (tested up to 50,000 transcripts)
  - Memory-efficient implementation

### Documentation

* Comprehensive function documentation with examples
* Introductory vignette with workflow demonstration
* Detailed README with troubleshooting guide

### Testing

* Unit tests covering all major functions
* Validation against synthetic data with known ground truth
* Reproducibility tests ensuring consistent results

### Known Limitations

* Requires at least one "Total" sample per condition-replicate combination
* Minimum of 10 transcripts recommended for reliable regression
* Statistical tests require at least 2 replicates per condition

### Future Development

We welcome contributions and feedback. Please report issues at:
https://github.com/Arnaroo/FracFixR/issues