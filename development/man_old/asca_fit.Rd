% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asca_fit.R
\name{asca_fit}
\alias{asca_fit}
\title{ASCA Fitting Workhorse Function}
\usage{
asca_fit(
  formula,
  data,
  subset,
  weights,
  na.action,
  family,
  permute = FALSE,
  perm.type = c("approximate", "exact"),
  unrestricted = FALSE,
  add_error = FALSE,
  aug_error = "denominator",
  use_ED = FALSE,
  pca.in = FALSE,
  coding = c("sum", "weighted", "reference", "treatment"),
  SStype = "II",
  REML = NULL
)
}
\arguments{
\item{formula}{Model formula accepting a single response (block) and predictors. See Details for more information.}

\item{data}{The data set to analyse.}

\item{subset}{Expression for subsetting the data before modelling.}

\item{weights}{Optional object weights.}

\item{na.action}{How to handle NAs (no action implemented).}

\item{family}{Error distributions and link function for Generalized Linear Models.}

\item{permute}{Perform approximate permutation testing, default = FALSE (numeric or TRUE = 1000).}

\item{perm.type}{Type of permutation: "approximate" (default) or "exact".}

\item{unrestricted}{Use unrestricted ANOVA decomposition (default = FALSE).}

\item{add_error}{Add error to LS means, e.g., for APCA.}

\item{aug_error}{Augment score matrices in backprojection. Default = "denominator"
(of F test), "residual" (force error term), nueric value (alpha-value in LiMM-PCA).}

\item{use_ED}{Use "effective dimensions" for score rescaling in LiMM-PCA.}

\item{pca.in}{Compress response before ASCA (number of components).}

\item{coding}{Effect coding: "sum" (default = sum-coding), "weighted", "reference", "treatment".}

\item{SStype}{Type of sum-of-squares: "I" = sequential, "II" (default) = last term, obeying marginality,
"III" = last term, not obeying marginality.}

\item{REML}{Parameter to mixlm: NULL (default) = sum-of-squares, TRUE = REML, FALSE = ML.}
}
\description{
This function is called by all ASCA related methods in this package. It is documented
so that one can have access to a richer set of parameters from the various methods or call this
function directly. The latter should be done with care as there are many possibilities and not all
have been used in publications or tested thoroughly.
}
