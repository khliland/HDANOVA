# Extracting the Extended Model Frame from a Formula or Fit

This function attempts to apply
[`model.frame`](https://rdrr.io/r/stats/model.frame.html) and extend the
result with columns of interactions.

## Usage

``` r
extended.model.frame(formula, data, ..., sep = ".")
```

## Arguments

- formula:

  a model formula or terms object or an R object.

- data:

  a data.frame, list or environment (see
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html)).

- ...:

  further arguments to pass to
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html).

- sep:

  separator in contraction of names for interactions (default = ".").

## Value

A [`data.frame`](https://rdrr.io/r/base/data.frame.html) that includes
everything a [`model.frame`](https://rdrr.io/r/stats/model.frame.html)
does plus interaction terms.

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
and
[`pcanova_plots`](https://khliland.github.io/HDANOVA/reference/pcanova_plots.md)

## Examples

``` r
dat <- data.frame(Y = c(1,2,3,4,5,6),
                  X = factor(LETTERS[c(1,1,2,2,3,3)]),
                  W = factor(letters[c(1,2,1,2,1,2)]))
extended.model.frame(Y ~ X*W, dat)
#>   Y X W X:W
#> 1 1 A a A.a
#> 2 2 A b A.b
#> 3 3 B a B.a
#> 4 4 B b B.b
#> 5 5 C a C.a
#> 6 6 C b C.b
```
