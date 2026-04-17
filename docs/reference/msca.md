# Multilevel Simultaneous Component Analysis - MSCA

This MSCA implementation assumes a single factor to be used as
between-individuals factor.

## Usage

``` r
msca(
  formula,
  data,
  contrasts = "contr.sum",
  permute = FALSE,
  perm.type = c("approximate", "exact"),
  ...
)
```

## Arguments

- formula:

  Model formula accepting a single response (block) and predictors. See
  Details for more information.

- data:

  The data set to analyse.

- contrasts:

  Effect coding: "sum" (default = sum-coding), "weighted", "reference",
  "treatment".

- permute:

  Number of permutations to perform (default = 1000).

- perm.type:

  Type of permutation to perform, either "approximate" or "exact"
  (default = "approximate").

- ...:

  Additional arguments to
  [`hdanova`](https://khliland.github.io/HDANOVA/reference/hdanova.md).

## Value

An `asca` object containing loadings, scores, explained variances, etc.
The object has associated plotting
([`asca_plots`](https://khliland.github.io/HDANOVA/reference/asca_plots.md))
and result
([`asca_results`](https://khliland.github.io/HDANOVA/reference/asca_results.md))
functions.

## References

- Smilde, A., Jansen, J., Hoefsloot, H., Lamers,R., Van Der Greef, J.,
  and Timmerman, M.(2005). ANOVA-Simultaneous Component Analysis (ASCA):
  A new tool for analyzing designed metabolomics data. Bioinformatics,
  21(13), 3043–3048.

- Liland, K.H., Smilde, A., Marini, F., and Næs,T. (2018). Confidence
  ellipsoids for ASCA models based on multivariate regression theory.
  Journal of Chemometrics, 32(e2990), 1–13.

## See also

Main methods:
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.md),
[`apca`](https://khliland.github.io/HDANOVA/reference/apca.md),
[`limmpca`](https://khliland.github.io/HDANOVA/reference/limmpca.md),
`msca`,
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
# Load candies data
data(candies)

# Basic MSCA model with a single factor
mod <- msca(assessment ~ candy, data=candies)
print(mod)
#> Multilevel Simultaneous Component Analysis fitted using 'lm' (Linear Model)
#> Call:
#> msca(formula = assessment ~ candy, data = candies)
summary(mod)
#> Multilevel Simultaneous Component Analysis fitted using 'lm' (Linear Model) 
#> - SS type II, sum coding, restricted model, least squares estimation, SSQ method: qr_regression 
#>          Sum.Sq. Expl.var.(%)
#> Between 33416.66        74.48
#> Within  11450.62        25.52

# Result plotting for first factor
loadingplot(mod, scatter=TRUE, labels="names")

scoreplot(mod)


# Within scores
scoreplot(mod, factor="within")


# Within scores per factor level
par.old <- par(mfrow=c(3,2), mar=c(4,4,2,1), mgp=c(2,0.7,0))
for(i in 1:length(mod$scores.within))
  scoreplot(mod, factor="within", within_level=i,
            main=paste0("Level: ", names(mod$scores.within)[i]),
            panel.first=abline(v=0,h=0,col="gray",lty=2))
par(par.old)


# Permutation testing
mod.perm <- asca(assessment ~ candy * assessor, data=candies, permute=TRUE)
summary(mod.perm)
#> Anova Simultaneous Component Analysis fitted using 'lm' (Linear Model) 
#> - SS type II, sum coding, restricted model, least squares estimation, SSQ method: qr_regression, 1000 permutations 
#>                 Sum.Sq. Expl.var.(%) p-value
#> candy          33416.66        74.48       0
#> assessor        1961.37         4.37       0
#> candy:assessor  3445.73         7.68       0
#> Residuals       6043.52        13.47      NA
```
