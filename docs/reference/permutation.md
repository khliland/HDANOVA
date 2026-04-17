# Permutation for SS-Type-Aligned HDANOVA

Permutation testing for
[`hdanova()`](https://khliland.github.io/HDANOVA/reference/hdanova.md)
objects using SS-type-aligned QR effect matrices. This function is
intended as an opt-in alternative to the legacy regression-based
`permutation()` workflow.

## Usage

``` r
permutation(
  object,
  permute = 1000,
  perm.type = c("approximate", "exact"),
  respect_SStype = NULL,
  unique.digits = 12,
  unique.frac = 0.95,
  exhaustive.warn = TRUE
)
```

## Arguments

- object:

  A `hdanova` object.

- permute:

  Number of permutations to perform (default = 1000).

- perm.type:

  Type of permutation to perform, either `"approximate"` or `"exact"`
  (default = `"approximate"`).

- respect_SStype:

  Logical or `NULL`. If `NULL` (default), use the `hdanova` object
  setting (`object$more$respect_SStype`). If `FALSE`, follow the legacy
  regression-based permutation logic. If `TRUE`, use SS-type-aligned QR
  contrasts. For REML/ML mixed models, this may yield permutation
  statistics that differ from fitted-model REML SSQ decompositions
  selected by `REML_ssq_method`.

- unique.digits:

  Number of digits used when rounding permutation SSQ values before
  checking uniqueness (default = 12). Set to `NULL` to disable this
  warning.

- unique.frac:

  Minimum fraction of unique rounded SSQ values required to avoid
  warning (default = 0.95). Set to `NULL` to disable this warning.

- exhaustive.warn:

  Logical; if `TRUE` (default), warn when exact permutation uses
  exhaustive enumeration with fewer permutations than requested.

## Value

An updated `hdanova` object with `permute` results.

## Details

The permutation statistics are computed from SS-type-aligned QR
reduced/full model contrasts rather than the legacy regression-based LS
matrices. Fixed models and mixed MoM models are supported. Approximate
permutation uses a relaxed global shuffle of observations; exact
permutation uses permissible block-restricted shuffles. For REML/ML
mixed models with `respect_SStype = TRUE`, a warning is issued to
highlight that permutation statistics and fitted-model REML SSQ
decompositions are based on different computational definitions.
