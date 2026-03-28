# Permutation Based MANOVA - PERMANOVA

Wrapper for the
[`adonis2`](https://vegandevs.github.io/vegan/reference/adonis.html)
function to allow ordinary formula input.

## Usage

``` r
permanova(formula, data, ...)
```

## Arguments

- formula:

  Model formula accepting a single response matrix and predictors. See
  details in
  [`adonis2`](https://vegandevs.github.io/vegan/reference/adonis.html).

- data:

  The data set to analyse.

- ...:

  Additional arguments to
  [`adonis2`](https://vegandevs.github.io/vegan/reference/adonis.html).

## Value

An ANOVA table with permutation-based p-values.

## Examples

``` r
data(caldana)
(pr <- permanova(compounds ~ light * time, caldana))
#> Permutation test for adonis under reduced model
#> Permutation: free
#> Number of permutations: 999
#> 
#> vegan::adonis2(formula = formula, data = data)
#>           Df SumOfSqs      R2      F Pr(>F)    
#> Model     27  0.79189 0.35355 2.2687  0.001 ***
#> Residual 112  1.44793 0.64645                  
#> Total    139  2.23982 1.00000                  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
