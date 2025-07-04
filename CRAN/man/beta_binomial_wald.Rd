% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fracfixr.r
\name{beta_binomial_wald}
\alias{beta_binomial_wald}
\title{beta_binomial_wald: Beta-Binomial Wald Test}
\usage{
beta_binomial_wald(counts, successes, annotation)
}
\arguments{
\item{counts}{Total count matrix}

\item{successes}{Success count matrix (not proportions)}

\item{annotation}{Sample metadata}
}
\value{
Data frame with test results for all transcripts
}
\description{
Implements Wald test using beta-binomial distribution to account for
overdispersion in count data. Useful when variance exceeds that expected
under binomial distribution.
}
\keyword{internal}
