% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{prc}
\alias{prc}
\title{Principal Response Curves}
\usage{
prc(formula, data, ...)
}
\arguments{
\item{formula}{Model formula accepting a single response (block) and predictors.
If no predictor is called 'time', time is assumed to be the second predictor.}

\item{data}{The data set to analyse.}

\item{...}{Additional arguments to \code{\link[vegan]{prc}}.}
}
\value{
An object of class \code{prc}.
}
\description{
Wrapper for the \code{\link[vegan]{prc}} function to allow for formula input.
}
\examples{
data(caldana)
(pr <- prc(compounds ~ light * time, caldana))
summary(pr)
}
\seealso{
Main methods: \code{\link{asca}}, \code{\link{apca}}, \code{\link{limmpca}}, \code{\link{msca}}, \code{\link{pcanova}}, \code{\link{prc}} and \code{\link{permanova}}.
Workhorse function underpinning most methods: \code{\link{hdanova}}.
Extraction of results and plotting: \code{\link{asca_results}}, \code{\link{asca_plots}}, \code{\link{pcanova_results}} and \code{\link{pcanova_plots}}
}
