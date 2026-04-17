# ANOVA Principal Component Analysis - APCA

APCA function for fitting ANOVA Principal Component Analysis models.

## Usage

``` r
apca(
  formula,
  data,
  add_error = TRUE,
  contrasts = "contr.sum",
  permute = FALSE,
  perm.type = c("approximate", "exact"),
  ...
)
```

## Arguments

- formula:

  Model formula accepting a single response (block) and predictors.

- data:

  The data set to analyse.

- add_error:

  Add error to LS means (default = TRUE).

- contrasts:

  Effect coding: "sum" (default = sum-coding), "weighted", "reference",
  "treatment".

- permute:

  Number of permutations to perform (default = 1000).

- perm.type:

  Type of permutation to perform, either "approximate" or "exact"
  (default = "approximate").

- ...:

  Additional parameters for the `hdanova` function.

## Value

An object of class `apca`, inheriting from the general `asca` class.
Further arguments and plots can be found in the
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.md)
documentation.

## References

Harrington, P.d.B., Vieira, N.E., Espinoza, J., Nien, J.K., Romero, R.,
and Yergey, A.L. (2005) Analysis of variance–principal component
analysis: A soft tool for proteomic discovery. Analytica chimica acta,
544 (1-2), 118–127.

## See also

Main methods:
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.md), `apca`,
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
data(candies)
ap <- apca(assessment ~ candy, data=candies)
scoreplot(ap)


# Numeric effects
candies$num <- eff <- 1:165
mod <- apca(assessment ~ candy + assessor + num, data=candies)
summary(mod)
#> Anova Principal Component Analysis fitted using 'lm' (Linear Model) 
#> - SS type II, sum coding, restricted model, least squares estimation, SSQ method: qr_regression 
#>            Sum.Sq. Expl.var.(%)
#> candy     32438.46        72.30
#> assessor   1823.59         4.06
#> num         101.17         0.23
#> Residuals  9388.08        20.92
scoreplot(mod, factor=3, gr.col=rgb(eff/max(eff), 1-eff/max(eff),0), pch.scores="x")
```
