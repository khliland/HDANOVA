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

## Examples

``` r
# Load candies data
data(candies)

# Basic HDANOVA model with two factors
mod <- hdanova(assessment ~ candy + assessor, data=candies)
mod <- sca(mod)
scoreplot(mod)

```
