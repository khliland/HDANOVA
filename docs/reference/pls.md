# Partial Least Squares (PLS) for HDANOVA

This function performs Partial Least Squares (PLS) on a `hdanova`.

## Usage

``` r
pls(object, ...)

# Default S3 method
pls(object, ...)
```

## Arguments

- object:

  A `hdanova` object.

- ...:

  Additional arguments (not used).

## Value

An updated `hdanova` object with PLS results.

## Details

For residuals, PCA is performed instead of PLS as there is no natural
response.

## See also

Main wrapper:
[`apls`](https://khliland.github.io/HDANOVA/reference/apls.md). Related
decomposition:
[`sca`](https://khliland.github.io/HDANOVA/reference/sca.md). Plotting
and summaries:
[`asca_plots`](https://khliland.github.io/HDANOVA/reference/asca_plots.md)
and
[`asca_results`](https://khliland.github.io/HDANOVA/reference/asca_results.md).

## Examples

``` r
# Load candies data
data(candies)

# Basic HDANOVA model with two factors
mod <- hdanova(assessment ~ candy + assessor, data=candies)
mod <- pls(mod)
scoreplot(mod)

```
