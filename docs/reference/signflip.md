# Flip signs of a component/factor combination in a SCA/PCA object

This function flips the sign of a selected component in a selected
factor of an `asca` object. This affects both scores, loadings and
projected data.

## Usage

``` r
signflip(object, factor, comp)
```

## Arguments

- object:

  `asca` object.

- factor:

  `integer/character` for selecting a model factor.

- comp:

  `integer` for selected component.

## Value

An `asca` object with the sign of the selected component flipped.

## Examples

``` r
# Load candies data
data(candies)

# Basic HDANOVA model with two factors
mod <- hdanova(assessment ~ candy + assessor, data=candies)
mod <- sca(mod)
old.par <- par(mfrow=c(1,2), mar=c(4,4,1,1))
scoreplot(mod, factor="candy")
loadingplot(mod, factor="candy")

par(old.par)

# Flip the sign of the first component of the candy factor
mod <- signflip(mod, factor="candy", comp=1)
old.par <- par(mfrow=c(1,2), mar=c(4,4,1,1))
scoreplot(mod, factor="candy")
loadingplot(mod, factor="candy")

par(old.par)
```
