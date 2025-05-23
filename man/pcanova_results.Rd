% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcanova_results.R
\name{pcanova_results}
\alias{pcanova_results}
\alias{summary.pcanova}
\alias{projections.pcanova}
\alias{print.summary.pcanova}
\alias{print.pcanova}
\title{PC-ANOVA Result Methods}
\usage{
\method{summary}{pcanova}(object, ...)

\method{print}{summary.pcanova}(x, digits = 2, ...)

\method{print}{pcanova}(x, ...)

\method{summary}{pcanova}(object, ...)
}
\arguments{
\item{object}{\code{pcanova} object.}

\item{...}{additional arguments to underlying methods.}

\item{x}{\code{pcanova} object.}

\item{digits}{\code{integer} number of digits for printing.}
}
\value{
Returns depend on method used, e.g. \code{projections.pcanova} returns projected samples,
\code{scores.pcanova} return scores, while print and summary methods return the object invisibly.
}
\description{
Standard result computation and extraction functions for ASCA (\code{\link{pcanova}}).
}
\details{
Usage of the functions are shown using generics in the examples in \code{\link{pcanova}}.
Explained variances are available (block-wise and global) through \code{blockexpl} and \code{print.rosaexpl}.
Object printing and summary are available through:
\code{print.pcanova} and \code{summary.pcanova}.
Scores and loadings have their own extensions of \code{scores()} and \code{loadings()} through
\code{scores.pcanova} and \code{loadings.pcanova}. Special to ASCA is that scores are on a
factor level basis, while back-projected samples have their own function in \code{projections.pcanova}.
}
\references{
Luciano G, Næs T. Interpreting sensory data by combining principal
component analysis and analysis of variance. Food Qual Prefer. 2009;20(3):167-175.
}
\seealso{
Main methods: \code{\link{asca}}, \code{\link{apca}}, \code{\link{limmpca}}, \code{\link{msca}}, \code{\link{pcanova}}, \code{\link{prc}} and \code{\link{permanova}}.
Workhorse function underpinning most methods: \code{\link{hdanova}}.
Extraction of results and plotting: \code{\link{asca_results}}, \code{\link{asca_plots}}, \code{\link{pcanova_results}} and \code{\link{pcanova_plots}}
}
