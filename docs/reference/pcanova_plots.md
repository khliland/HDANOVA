# PC-ANOVA Result Methods

Various plotting procedures for
[`pcanova`](https://khliland.github.io/HDANOVA/reference/pcanova.md)
objects.

## Usage

``` r
# S3 method for class 'pcanova'
scoreplot(object, factor = 1, comps = 1:2, col = "factor", ...)
```

## Arguments

- object:

  `pcanova` object.

- factor:

  `integer/character` for selecting a model factor.

- comps:

  `integer` vector of selected components.

- col:

  `character` for selecting a factor to use for colouring (default =
  first factor) or ordinary colour specifications.

- ...:

  additional arguments to underlying methods.

## Value

The plotting routines have no return.

## Details

Usage of the functions are shown using generics in the examples in
[`pcanova`](https://khliland.github.io/HDANOVA/reference/pcanova.md).
Plot routines are available as `scoreplot.pcanova` and
`loadingplot.pcanova`.

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
[`pcanova_results`](https://khliland.github.io/HDANOVA/reference/pcanova_results.md)
and `pcanova_plots`
