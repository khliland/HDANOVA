# ASCA Fitting Workhorse Function

This function is called by all ASCA related methods in this package. It
is documented so that one can have access to a richer set of parameters
from the various methods or call this function directly. The latter
should be done with care as there are many possibilities and not all
have been used in publications or tested thoroughly.

## Usage

``` r
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
  contrasts = "contr.sum",
  coding,
  equal_baseline = FALSE,
  SStype = "II",
  REML = NULL
)
```

## Arguments

- formula:

  Model formula accepting a single response (block) and predictors. See
  Details for more information.

- data:

  The data set to analyse.

- subset:

  Expression for subsetting the data before modelling.

- weights:

  Optional object weights.

- na.action:

  How to handle NAs (no action implemented).

- family:

  Error distributions and link function for Generalized Linear Models.

- permute:

  Perform approximate permutation testing, default = FALSE (numeric or
  TRUE = 1000 permutations).

- perm.type:

  Type of permutation: "approximate" (default) or "exact".

- unrestricted:

  Use unrestricted ANOVA decomposition (default = FALSE).

- add_error:

  Add error to LS means, e.g., for APCA.

- aug_error:

  Augment score matrices in backprojection. Default = "denominator" (of
  F test), "residual" (force error term), nueric value (alpha-value in
  LiMM-PCA).

- use_ED:

  Use "effective dimensions" for score rescaling in LiMM-PCA.

- pca.in:

  Compress response before ASCA (number of components).

- contrasts:

  Effect coding: "sum" (default = sum-coding), "weighted", "reference",
  "treatment".

- coding:

  Defunct. Use 'contrasts' instead.

- equal_baseline:

  Experimental: Set to `TRUE` to let interactions, where a main effect
  is missing, e.g., a nested model, be handled with the same baseline as
  a cross effect model. If `TRUE` the corresponding interactions will be
  put in quotation marks and included in the `model.frame`.

- SStype:

  Type of sum-of-squares: "I" = sequential, "II" (default) = last term,
  obeying marginality, "III" = last term, not obeying marginality.

- REML:

  Parameter to mixlm: NULL (default) = sum-of-squares, TRUE = REML,
  FALSE = ML.

## Value

An `asca` object containing loadings, scores, explained variances, etc.
The object has associated plotting
([`asca_plots`](https://khliland.github.io/HDANOVA/reference/asca_plots.md))
and result
([`asca_results`](https://khliland.github.io/HDANOVA/reference/asca_results.md))
functions.
