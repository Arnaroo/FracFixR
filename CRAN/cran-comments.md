## Test environments

* local Ubuntu 22.04, R 4.3.2
* win-builder (devel and release)
* macOS (via GitHub Actions), R release
* Windows (via GitHub Actions), R release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Initial submission

This is the first submission of FracFixR to CRAN.

## Package purpose

FracFixR provides a novel statistical framework for analyzing RNA-seq data from fractionation experiments. It addresses a significant gap in current bioinformatics tools by properly accounting for the compositional nature of fractionated samples and estimating unrecoverable material ("lost" fraction).

## Dependencies

All dependencies are well-established packages available on CRAN or Bioconductor. We have kept dependencies minimal while ensuring full functionality.

## Examples

All examples run in < 5 seconds on the test environments. Examples using parallel processing are wrapped in \dontrun{} where appropriate.

## License

The package is released under CC BY-NC-ND 4.0 license as noted in the DESCRIPTION file. This is appropriate for academic software with specific usage restrictions.

## Bioconductor packages

The package depends on EnhancedVolcano from Bioconductor. This is properly declared in the DESCRIPTION file.

## Test coverage

The package includes comprehensive unit tests achieving >90% code coverage for critical functions.