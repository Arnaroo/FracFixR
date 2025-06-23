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