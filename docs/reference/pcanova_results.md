# PC-ANOVA Result Methods

Standard result computation and extraction functions for ASCA
([`pcanova`](https://khliland.github.io/HDANOVA/reference/pcanova.md)).

## Usage

``` r
# S3 method for class 'pcanova'
summary(object, ...)

# S3 method for class 'summary.pcanova'
print(x, digits = 2, ...)

# S3 method for class 'pcanova'
print(x, ...)

# S3 method for class 'pcanova'
summary(object, ...)
```

## Arguments

- object:

  `pcanova` object.

- ...:

  additional arguments to underlying methods.

- x:

  `pcanova` object.

- digits:

  `integer` number of digits for printing.

## Value

Returns depend on method used, e.g. `projections.pcanova` returns
projected samples, `scores.pcanova` return scores, while print and
summary methods return the object invisibly.

## Details

Usage of the functions are shown using generics in the examples in
[`pcanova`](https://khliland.github.io/HDANOVA/reference/pcanova.md).
Explained variances are available (block-wise and global) through
`blockexpl` and `print.rosaexpl`. Object printing and summary are
available through: `print.pcanova` and `summary.pcanova`. Scores and
loadings have their own extensions of
[`scores()`](https://khliland.github.io/pls/reference/scores.html) and
[`loadings()`](https://khliland.github.io/pls/reference/scores.html)
through `scores.pcanova` and `loadings.pcanova`. Special to ASCA is that
scores are on a factor level basis, while back-projected samples have
their own function in `projections.pcanova`.

## References

Luciano G, Næs T. Interpreting sensory data by combining principal
component analysis and analysis of variance. Food Qual Prefer.
2009;20(3):167-175.

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
Extraction of results and plotting:
[`asca_results`](https://khliland.github.io/HDANOVA/reference/asca_results.md),
[`asca_plots`](https://khliland.github.io/HDANOVA/reference/asca_plots.md),
`pcanova_results` and
[`pcanova_plots`](https://khliland.github.io/HDANOVA/reference/pcanova_plots.md)
