% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fracfixr.r
\name{PlotComparison}
\alias{PlotComparison}
\title{PlotComparison: Create Volcano Plot for Differential Results}
\usage{
PlotComparison(DiffPropResult, Conditions = NULL, Types = NULL, cutoff = NULL)
}
\arguments{
\item{DiffPropResult}{Output from DiffPropTest() function}

\item{Conditions}{Character vector of conditions being compared}

\item{Types}{Character vector of fraction types analyzed}

\item{cutoff}{Optional y-axis maximum for plot}
}
\value{
EnhancedVolcano plot object
}
\description{
Generates an enhanced volcano plot showing transcripts with significant
differential proportions between conditions.
}
\examples{
\dontrun{
# Create volcano plot
volcano <- PlotComparison(diff_results, 
                         Conditions = c("Control", "Treatment"),
                         Types = "Heavy_Polysome")
}

}
