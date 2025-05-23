% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signflip.R
\name{signflip}
\alias{signflip}
\title{Flip signs of a component/factor combination in a SCA/PCA object}
\usage{
signflip(object, factor, comp)
}
\arguments{
\item{object}{\code{asca} object.}

\item{factor}{\code{integer/character} for selecting a model factor.}

\item{comp}{\code{integer} for selected component.}
}
\value{
An \code{asca} object with the sign of the selected component flipped.
}
\description{
This function flips the sign of a selected component in a selected factor of an \code{asca} object.
This affects both scores, loadings and projected data.
}
\examples{
# Load candies data
data(candies)

# Basic HDANOVA model with two factors
mod <- hdanova(assessment ~ candy + assessor, data=candies)
mod <- sca(mod)
old.par <- par(mfrow=c(1,2), mar=c(4,4,1,1))
scoreplot(mod, factor="candy")
loadingplot(mod, factor="candy")
par(old.par)

# Flip the sign of the first component of the candy factor
mod <- signflip(mod, factor="candy", comp=1)
old.par <- par(mfrow=c(1,2), mar=c(4,4,1,1))
scoreplot(mod, factor="candy")
loadingplot(mod, factor="candy")
par(old.par)
}
