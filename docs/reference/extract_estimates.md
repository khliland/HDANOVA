# Extract estimates for a given factor combination

Extracts and sums the LS estimates for a given factor combination from
an object of class `hdanova`. If `add_residuals` is `TRUE`, the
residuals are added to the LS estimates. If `substract` is `TRUE`, the
returned matrix is the data with chosen estimates subtracted.

## Usage

``` r
extract_estimates(object, factors, subtract = FALSE, add_residuals = FALSE)
```

## Arguments

- object:

  `asca` object.

- factors:

  `vector` of factor names or numbers.

- subtract:

  `logical` subtract the estimates from the data (default = FALSE).

- add_residuals:

  `logical` add residuals to the estimates (default = FALSE).

## Value

A matrix of the extracted estimates.

## Examples

``` r
# Load candies data
data(candies)

# Basic HDANOVA model with two factors and interaction
mod <- hdanova(assessment ~ candy * assessor, data=candies)

# Extract estimates for the interaction
inter <- extract_estimates(mod, c("candy:assessor"))

# Visualize the interaction effect
image(t(inter), main="Interaction effect", xlab="Attribute", ylab="Sample")

```
