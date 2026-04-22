# Simultaneous Component Analysis

This function performs Simultaneous Component Analysis (SCA) on a
`hdanova` object.

## Usage

``` r
sca(object)
```

## Arguments

- object:

  A `hdanova` object.

## Value

An updated `hdanova` object with SCA results.

## See also

Model constructors and wrappers:
[`hdanova`](https://khliland.github.io/HDANOVA/reference/hdanova.md),
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.md),
[`apca`](https://khliland.github.io/HDANOVA/reference/apca.md),
[`limmpca`](https://khliland.github.io/HDANOVA/reference/limmpca.md) and
[`msca`](https://khliland.github.io/HDANOVA/reference/msca.md). Plotting
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
mod <- sca(mod)
scoreplot(mod)

```
