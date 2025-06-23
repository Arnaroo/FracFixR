# Contributing to FracFixR

Thank you for your interest in contributing to FracFixR! As scientific software, we maintain strict standards to ensure reproducibility, accuracy, and reliability of all analyses. Please read these guidelines carefully before submitting contributions.

## Table of Contents

- [License Agreement](#license-agreement)
- [Code of Conduct](#code-of-conduct)
- [Scientific Integrity](#scientific-integrity)
- [How to Contribute](#how-to-contribute)
- [Development Setup](#development-setup)
- [Contribution Standards](#contribution-standards)
- [Testing Requirements](#testing-requirements)
- [Documentation Standards](#documentation-standards)
- [Submission Process](#submission-process)
- [Review Process](#review-process)

## License Agreement

FracFixR is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License (CC BY-NC-ND 4.0).

**Important**: By contributing to FracFixR, you agree that:

1. Your contributions will be licensed under the same CC BY-NC-ND 4.0 license
2. You have the right to submit the contribution
3. Your contribution does not infringe on any third-party rights
4. You understand that your contribution may be included in academic publications

**Note on Derivatives**: While the CC BY-NC-ND license restricts derivative works, we consider contributions that maintain the core scientific methodology and improve the implementation as enhancements rather than derivatives. However, any fundamental changes to the statistical framework require explicit permission from the original authors.

## Code of Conduct

We are committed to fostering an open, welcoming, and harassment-free environment. All contributors must:

1. Use welcoming and inclusive language
2. Respect differing viewpoints and experiences
3. Accept constructive criticism gracefully
4. Focus on what's best for the scientific community
5. Show empathy towards other community members

Unacceptable behavior includes harassment, discriminatory language, and other conduct that creates an unsafe environment. Violations may result in removal from the project.

## Scientific Integrity

As scientific software, FracFixR must maintain the highest standards of accuracy and reproducibility:

### Core Principles

1. **Reproducibility First**: All code must produce identical results given the same input and random seed
2. **Mathematical Accuracy**: Any changes to statistical methods must be mathematically sound and peer-reviewed
3. **Transparency**: All algorithms and methods must be clearly documented
4. **Validation**: New features must be validated against known results or published data

### Protected Components

The following components are core to the scientific methodology and require special approval to modify:

- NNLS regression implementation in `ProcessReplicate()`
- Transcript selection criteria (70-96% quantile range)
- Statistical test implementations in `DiffPropTest()`
- Proportion calculation formulas

## How to Contribute

### Types of Contributions Welcome

1. **Bug Fixes**: Corrections to ensure accurate results
2. **Performance Improvements**: Optimizations that maintain identical outputs
3. **Documentation**: Clarifications, examples, and tutorials
4. **Visualization Enhancements**: Improved plots that better communicate results
5. **Test Cases**: Additional validation scenarios
6. **Platform Compatibility**: Ensuring consistent behavior across systems

### Contributions Requiring Special Review

1. **New Statistical Methods**: Must include mathematical justification and citations
2. **Algorithm Changes**: Must demonstrate equivalence or improvement with validation
3. **Default Parameter Changes**: Must be justified with empirical evidence
4. **Core Function Modifications**: Require extensive testing and validation

### Contributions Not Accepted

1. Changes that alter the fundamental statistical framework without peer review
2. Features that compromise reproducibility
3. Non-standard dependencies that affect portability
4. Modifications that break backward compatibility without justification

## Development Setup

### Prerequisites

```bash
# R version ≥ 4.0.0
R --version

# Install development dependencies
R -e "install.packages(c('devtools', 'testthat', 'covr', 'lintr', 'roxygen2'))"
```

### Fork and Clone

```bash
# Fork the repository on GitHub, then:
git clone https://github.com/YOUR_USERNAME/FracFixR.git
cd FracFixR
git remote add upstream https://github.com/Arnaroo/FracFixR.git
```

### Development Environment

```r
# Load development environment
library(devtools)
load_all()

# Run tests
test()

# Check package
check()
```

## Contribution Standards

### Code Style

We follow the tidyverse style guide with specific requirements for scientific computing:

```r
# GOOD: Clear variable names with units
transcript_counts <- matrix(...)
fraction_proportion <- 0.7

# BAD: Ambiguous names
x <- matrix(...)
p <- 0.7

# GOOD: Explicit parameter validation
stopifnot(is.matrix(counts),
          is.numeric(counts),
          all(counts >= 0))

# GOOD: Documented magic numbers
s1 <- quantile(TotalSum, 0.7)  # Lower bound: 70th percentile
s2 <- quantile(TotalSum, 0.96) # Upper bound: 96th percentile
```

### Reproducibility Requirements

1. **Set Seeds**: Always set random seeds for any stochastic process

```r
set.seed(42)  # Use consistent seeds in examples
```

2. **Floating Point Comparisons**: Use appropriate tolerance

```r
# GOOD
all.equal(result, expected, tolerance = 1e-10)

# BAD
result == expected
```

3. **Platform Independence**: Avoid platform-specific code

```r
# GOOD
file.path("data", "counts.csv")

# BAD
"data/counts.csv"  # Assumes Unix-style paths
```

### Function Requirements

All functions must include:

```r
#' Function Title (Active Voice)
#'
#' @description
#' Detailed description of what the function does, including mathematical
#' basis if applicable.
#'
#' @param x Description with type and constraints
#' @param method Character string, one of "GLM", "Logit", "Wald"
#'
#' @return Description of return value with structure
#'
#' @details
#' Mathematical formulation:
#' \deqn{Y_{ij} = \alpha_0 + \sum_{f} \alpha_f X_{ijf} + \epsilon}
#'
#' @examples
#' # Reproducible example
#' set.seed(123)
#' data <- generate_test_data()
#' result <- function_name(data)
#'
#' @references
#' Cleynen et al. (2024) FracFixR: A compositional framework...
#'
#' @export
function_name <- function(x, method = c("GLM", "Logit", "Wald")) {
  # Input validation
  stopifnot(is.matrix(x), is.numeric(x))
  method <- match.arg(method)
  
  # Implementation
  ...
}
```

## Testing Requirements

### Test Coverage

- Minimum 90% code coverage for new functions
- 100% coverage for statistical computations

### Test Structure

```r
# tests/testthat/test-function_name.R
context("function_name")

test_that("function_name produces correct output", {
  # Setup
  set.seed(123)
  test_data <- generate_known_data()
  expected <- known_result()
  
  # Test
  result <- function_name(test_data)
  
  # Validate
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("function_name handles edge cases", {
  # Test with minimal data
  expect_error(function_name(matrix(1)), "at least 2")
  
  # Test with missing values
  data_na <- test_data
  data_na[1, 1] <- NA
  expect_error(function_name(data_na), "missing values")
})

test_that("function_name is reproducible", {
  set.seed(42)
  result1 <- function_name(test_data)
  
  set.seed(42)
  result2 <- function_name(test_data)
  
  expect_identical(result1, result2)
})
```

### Validation Requirements

For statistical methods, include:

1. Validation against published results
2. Simulation studies showing statistical properties
3. Comparison with alternative implementations

```r
# tests/testthat/test-validation.R
test_that("NNLS regression matches published results", {
  # Use data from Cleynen et al. 2024 supplementary
  published_data <- read.csv("tests/data/published_example.csv")
  published_result <- 0.847  # From paper
  
  our_result <- ProcessReplicate(published_data, ...)
  
  expect_equal(our_result$coefficient, published_result, tolerance = 0.001)
})
```

## Documentation Standards

### README Updates

- Add examples for new features
- Update function reference
- Include citations for new methods

### Vignette Requirements

For substantial features, create a vignette:

```r
# vignettes/new_feature.Rmd
---
title: "Using New Feature in FracFixR"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using New Feature}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```

### Mathematical Documentation

- Use LaTeX notation for formulas
- Include derivations for new methods
- Cite relevant literature

```r
#' The optimization problem:
#' \deqn{\min_{\alpha \geq 0} ||Y - X\alpha||_2^2}
#' 
#' Where \eqn{Y} is the total RNA vector and \eqn{X} is the fraction matrix.
```

## Submission Process

### 1. Create an Issue First

Before starting work:

```markdown
**Issue Title**: [BUG/FEATURE/DOC] Clear description

**Description**: 
- What problem does this solve?
- How does it maintain reproducibility?
- Any mathematical basis?

**Validation Plan**:
- How will you test this?
- What published data can validate it?
```

### 2. Branch Naming

```bash
# Bug fixes
git checkout -b fix/issue-number-description

# Features
git checkout -b feature/issue-number-description

# Documentation
git checkout -b docs/issue-number-description
```

### 3. Commit Standards

```bash
# Good commit messages
git commit -m "fix: correct NNLS convergence for edge case (#42)

- Add tolerance parameter to nnls() call
- Validate against published dataset
- Add test case for low-count transcripts"

# Include issue number
git commit -m "feat: add bootstrap confidence intervals (#38)"
```

### 4. Pull Request Template

```markdown
## Description
Brief description of changes

## Issue
Fixes #[issue number]

## Type of Change
- [ ] Bug fix (non-breaking change fixing an issue)
- [ ] New feature (non-breaking change adding functionality)
- [ ] Breaking change (fix or feature causing existing functionality to change)
- [ ] Documentation update

## Validation
- [ ] Passes all existing tests
- [ ] New tests added (coverage ≥90%)
- [ ] Validated against published results
- [ ] Reproducible (same seed = same results)

## Mathematical Basis (if applicable)
Description or citation for any mathematical changes

## Checklist
- [ ] Code follows project style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] No warnings in R CMD check
- [ ] Version number updated if needed
```

## Review Process

### Review Criteria

Pull requests are evaluated based on:

1. **Scientific Accuracy**: Do changes maintain or improve accuracy?
2. **Reproducibility**: Are results identical with same inputs?
3. **Code Quality**: Is code clear, documented, and efficient?
4. **Test Coverage**: Are all changes thoroughly tested?
5. **Documentation**: Is usage clear to end users?

### Review Timeline

- Initial response: Within 7 days
- Review completion: Within 14 days for simple changes
- Complex changes: May require additional review from domain experts

### Merge Requirements

1. All tests pass
2. No decrease in test coverage
3. R CMD check passes with no warnings
4. At least one maintainer approval
5. For statistical changes: Review by statistical expert

## Questions?

For questions about contributing:

1. Check existing issues and discussions
2. Read the documentation thoroughly
3. Open a discussion for general questions
4. Contact maintainers for specific guidance

Thank you for helping improve FracFixR and advancing reproducible science!

*These guidelines ensure FracFixR remains a reliable tool for the scientific community. We appreciate your understanding and commitment to these standards.*
