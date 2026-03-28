# Permutation for HDANOVA

Permutation testing for HDANOVA. This function performes permutation
testing for the effects in the HDANOVA model and adds them to the
`hdanova` object.

## Usage

``` r
permutation(
  object,
  permute = 1000,
  perm.type = c("approximate", "exact"),
  unique.digits = 12,
  unique.frac = 0.95,
  exhaustive.warn = TRUE
)
```

## Arguments

- object:

  A `hdanova` object.

- permute:

  Number of permutations to perform (default = 10000).

- perm.type:

  Type of permutation to perform, either "approximate" or "exact"
  (default = "approximate").

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

An updated `hdanova` object with permutation results.

## Details

The function supports both approximate and exact permutation testing.
Approximate testing randomly permutes the data and calculates the SSQ
for each permutation, while exact testing permutes the data according to
the exchangeable units defined by the model structure, ensuring that the
permutations respect the dependencies in the data. The current
implementation uses the regression model for estimation, meaning that
unbalanced data will affect the SSQ values and thus the permutation
distribution.

## Examples

``` r
# Load candies data
data(candies)

## Basic HDANOVA model with two factors
mod <- hdanova(assessment ~ candy + assessor, data=candies)

## Approximate permutation
modApprox <- permutation(mod)
summary(modApprox)
#> High-Dimensional Analysis of Variance fitted using 'lm' (Linear Model) 
#> - SS type II, sum coding, restricted model, least squares estimation, 1000 permutations 
#>            Sum.Sq. Expl.var.(%) p-value
#> candy     33416.66        74.48       0
#> assessor   1961.37         4.37       0
#> Residuals  9489.25        21.15      NA

# Plot permutation distribution for "candy" effect
permutationplot(modApprox, factor="candy")


## Exact permutation (warning if too few exchangeable units)
modExact  <- permutation(mod, perm.type="exact")
summary(modExact)
#> High-Dimensional Analysis of Variance fitted using 'lm' (Linear Model) 
#> - SS type II, sum coding, restricted model, least squares estimation, 1000 permutations 
#>            Sum.Sq. Expl.var.(%) p-value
#> candy     33416.66        74.48       0
#> assessor   1961.37         4.37       0
#> Residuals  9489.25        21.15      NA

# Reduced candy data (first two levels of each effect, two replicates)
reduced_candies <- candies[candies$candy %in% levels(candies$candy)[1:2] &
     candies$assessor %in% levels(candies$assessor)[1:2], ][-seq(1,12,by=3),]
mod_reduced <- hdanova(assessment ~ candy + assessor, data=reduced_candies)
mod_reducedApprox <- permutation(mod_reduced)
#> Warning: Effect 'candy' produced only 36 unique SSQ values (rounded to 12 decimals) out of 1,000 evaluated permutations.
#> Effect 'assessor' produced only 35 unique SSQ values (rounded to 12 decimals) out of 1,000 evaluated permutations.
mod_reducedExact  <- permutation(mod_reduced, perm.type="exact")
#> Warning: Effect 'candy' uses exhaustive permutation with 576 unique orders (requested 1,000).
#> Effect 'candy' produced only 19 unique SSQ values (rounded to 12 decimals) in exhaustive permutation.
#> Effect 'assessor' uses exhaustive permutation with 576 unique orders (requested 1,000).
#> Effect 'assessor' produced only 18 unique SSQ values (rounded to 12 decimals) in exhaustive permutation.
par.old <- par(mfrow=c(2,1))
permutationplot(mod_reducedApprox, factor="assessor", main="Approximate permutation")
permutationplot(mod_reducedExact, factor="assessor", main="Exact permutation")

par(par.old)

# Check how many exchangeable units were available (minimum per effect)
lapply(modExact$permute$exchangeable, function(x)min(x$block_sizes))
#> $candy
#> [1] 15
#> 
#> $assessor
#> [1] 33
#> 

# Caldana data (combined effects and exact permutation)
data(caldana)
mod.comb <- asca(compounds ~ time + comb(light + time:light), data=caldana)
mod.comb <- permutation(mod.comb, perm.type="exact")
summary(mod.comb)
#> Anova Simultaneous Component Analysis fitted using 'lm' (Linear Model) 
#> - SS type II, sum coding, restricted model, least squares estimation, 1000 permutations 
#>                  Sum.Sq. Expl.var.(%) p-value
#> time              154.58         9.69       0
#> light+time:light  349.64        21.92       0
#> Residuals        1091.14        68.39      NA
```
