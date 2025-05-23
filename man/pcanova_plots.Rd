% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcanova_plots.R
\name{pcanova_plots}
\alias{pcanova_plots}
\alias{scoreplot.pcanova}
\alias{loadingplot.pcanova}
\title{PC-ANOVA Result Methods}
\usage{
\method{scoreplot}{pcanova}(object, factor = 1, comps = 1:2, col = "factor", ...)
}
\arguments{
\item{object}{\code{pcanova} object.}

\item{factor}{\code{integer/character} for selecting a model factor.}

\item{comps}{\code{integer} vector of selected components.}

\item{col}{\code{character} for selecting a factor to use for colouring
(default = first factor) or ordinary colour specifications.}

\item{...}{additional arguments to underlying methods.}
}
\value{
The plotting routines have no return.
}
\description{
Various plotting procedures for \code{\link{pcanova}} objects.
}
\details{
Usage of the functions are shown using generics in the examples in \code{\link{pcanova}}.
Plot routines are available as
\code{scoreplot.pcanova} and \code{loadingplot.pcanova}.
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
