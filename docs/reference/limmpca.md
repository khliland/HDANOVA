# Linear Mixed Model PCA

This function mimics parts of the LiMM-PCA framework, combining ASCA+
and linear mixed models to analyse high-dimensional designed data. The
default is to use REML estimation and scaling of the backprojected
errors. See examples for alternatives.

## Usage

``` r
limmpca(
  formula,
  data,
  pca.in = 5,
  aug_error = 0.05,
  use_ED = FALSE,
  REML = TRUE,
  contrasts = "contr.sum",
  permute = FALSE,
  perm.type = c("approximate", "exact"),
  SStype = "III",
  ...
)
```

## Arguments

- formula:

  Model formula accepting a single response (block) and predictors. See
  Details for more information.

- data:

  The data set to analyse.

- pca.in:

  Compress response before ASCA (number of components), default = 5.

- aug_error:

  Error term of model ("denominator", "residual", numeric alpha-value).
  The latter implies the first with a scaling factor.

- use_ED:

  Use Effective Dimensions instead of degrees of freedom when scaling.

- REML:

  Use restricted maximum likelihood estimation. Alternatives: TRUE
  (default), FALSE (ML), NULL (least squares).

- contrasts:

  Effect coding: "sum" (default = sum-coding), "weighted", "reference",
  "treatment".

- permute:

  Number of permutations to perform (default = 1000).

- perm.type:

  Type of permutation to perform, either "approximate" or "exact"
  (default = "approximate").

- SStype:

  Type of sum-of-squares: "I" = sequential, "II" = last term, obeying
  marginality, "III" (default) = last term, not obeying marginality.

- ...:

  Additional arguments to
  [`hdanova`](https://khliland.github.io/HDANOVA/reference/hdanova.md).

## Value

An object of class `limmpca`, inheriting from the general `asca` class.

## Details

The Sum of Squares for the model is dependent on the SStype of the
model. For SStype = "I" and SStype = "II" the SSQ is based on LLR
(possibly inflating large contributions), while it is directly estimated
from the model for SStype = "III". SStype = "III" is the default for
LiMM-PCA and should be combined with sum coding. Sum of Squares for the
random effects are based on the variance components.

## References

- Martin, M. and Govaerts, B. (2020). LiMM-PCA: Combining ASCA+ and
  linear mixed models to analyse high-dimensional designed data. Journal
  of Chemometrics, 34(6), e3232.

## See also

Main methods:
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.md),
[`apca`](https://khliland.github.io/HDANOVA/reference/apca.md),
`limmpca`,
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
# Load candies data
data(candies)

# Default LiMM-PCA model with two factors and interaction, 5 PCA components
mod <- limmpca(assessment ~ candy*r(assessor), data=candies)
#> boundary (singular) fit: see help('isSingular')
summary(mod)
#> LiMM-PCA fitted using 'lmm' (Linear Mixed Model) 
#> - SS type III, sum coding, restricted model, REML estimation, SSQ method: exact_refit 
#>                 Sum.Sq. Expl.var.(%)
#> candy          33394.40        78.00
#> candy:assessor   621.64         1.45
#> assessor         793.50         1.85
#> Residuals       6162.61        14.39
scoreplot(mod, factor = "candy")


# LiMM-PCA with least squares estimation and 8 PCA components
modLS <- limmpca(assessment ~ candy*r(assessor), data=candies, REML=NULL, pca.in=8)
summary(modLS)
#> LiMM-PCA fitted using 'lmm' (Linear Mixed Model) 
#> - SS type III, sum coding, restricted model, least squares estimation, SSQ method: qr_regression 
#>                 Sum.Sq. Expl.var.(%)
#> candy          33415.98        74.73
#> assessor        1948.75         4.36
#> candy:assessor  3419.04         7.65
#> Residuals       5934.46        13.27
scoreplot(modLS, factor = "candy")


# Load Caldana data
data(caldana)

# Combining effects in LiMM-PCA (assuming light is a random factor)
mod.comb <- limmpca(compounds ~ time + comb(r(light) + r(time:light)), data=caldana, pca.in=8)
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
summary(mod.comb)
#> LiMM-PCA fitted using 'lmm' (Linear Mixed Model) 
#> - SS type III, sum coding, restricted model, REML estimation, SSQ method: exact_refit 
#>                  Sum.Sq. Expl.var.(%)
#> time              129.51        11.20
#> light+time:light   83.89         7.25
#> Residuals         833.20        72.05
```
