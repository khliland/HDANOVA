# ASCA Result Methods

Standard result computation and extraction functions for ASCA
([`asca`](https://khliland.github.io/HDANOVA/reference/asca.md)).

## Usage

``` r
# S3 method for class 'hdanova'
print(x, ...)

# S3 method for class 'hdanova'
summary(object, extended = TRUE, df = FALSE, ...)

# S3 method for class 'summary.hdanova'
print(x, digits = 2, ...)

# S3 method for class 'asca'
loadings(object, factor = 1, ...)

# S3 method for class 'asca'
scores(object, factor = 1, ...)

projections(object, ...)

# S3 method for class 'asca'
projections(object, factor = 1, ...)
```

## Arguments

- x:

  `asca` object.

- ...:

  additional arguments to underlying methods.

- object:

  `asca` object.

- extended:

  Extended output in summary (default = TRUE).

- df:

  Show degrees of freedom in summary (default = FALSE).

- digits:

  `integer` number of digits for printing.

- factor:

  `integer/character` for selecting a model factor.

## Value

Returns depend on method used, e.g. `projections.asca` returns projected
samples, `scores.asca` return scores, while print and summary methods
return the object invisibly.

## Details

Usage of the functions are shown using generics in the examples in
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.md).
Explained variances are available (block-wise and global) through
`blockexpl` and `print.rosaexpl`. Object printing and summary are
available through: `print.asca` and `summary.asca`. Scores and loadings
have their own extensions of
[`scores()`](https://khliland.github.io/pls/reference/scores.html) and
[`loadings()`](https://khliland.github.io/pls/reference/scores.html)
through `scores.asca` and `loadings.asca`. Special to ASCA is that
scores are on a factor level basis, while back-projected samples have
their own function in `projections.asca`.

## References

- Smilde, A., Jansen, J., Hoefsloot, H., Lamers,R., Van Der Greef, J.,
  and Timmerman, M.(2005). ANOVA-Simultaneous Component Analysis (ASCA):
  A new tool for analyzing designed metabolomics data. Bioinformatics,
  21(13), 3043–3048.

- Liland, K.H., Smilde, A., Marini, F., and Næs,T. (2018). Confidence
  ellipsoids for ASCA models based on multivariate regression theory.
  Journal of Chemometrics, 32(e2990), 1–13.

- Martin, M. and Govaerts, B. (2020). LiMM-PCA: Combining ASCA+ and
  linear mixed models to analyse high-dimensional designed data. Journal
  of Chemometrics, 34(6), e3232.

## See also

Main methods:
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.md),
[`apca`](https://khliland.github.io/HDANOVA/reference/apca.md),
[`limmpca`](https://khliland.github.io/HDANOVA/reference/limmpca.md),
[`msca`](https://khliland.github.io/HDANOVA/reference/msca.md),
[`pcanova`](https://khliland.github.io/HDANOVA/reference/pcanova.md),
[`prc`](https://khliland.github.io/HDANOVA/reference/prc.md) and
[`permanova`](https://khliland.github.io/HDANOVA/reference/permanova.md).
Workhorse function underpinning most methods:
[`hdanova`](https://khliland.github.io/HDANOVA/reference/hdanova.md).
Extraction of results and plotting: `asca_results`,
[`asca_plots`](https://khliland.github.io/HDANOVA/reference/asca_plots.md),
[`pcanova_results`](https://khliland.github.io/HDANOVA/reference/pcanova_results.md)
and
[`pcanova_plots`](https://khliland.github.io/HDANOVA/reference/pcanova_plots.md)

## Examples

``` r
# For end-to-end examples using summary(), scores(), loadings(), and
# projections(), see ?asca and ?apca.
```
