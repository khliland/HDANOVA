# Biplot for ASCA models

Biplot for ASCA models

## Usage

``` r
# S3 method for class 'asca'
biplot(
  x,
  factor = 1,
  comps = 1:2,
  xlim = NULL,
  ylim = NULL,
  col = "darkgray",
  expand = 1,
  labels,
  legendpos,
  ...
)
```

## Arguments

- x:

  `asca` object.

- factor:

  Factor number or name.

- comps:

  `integer` vector of selected components.

- xlim:

  `numeric` vector of length 2 for x-axis limits of the loadings.

- ylim:

  `numeric` vector of length 2 for y-axis limits of the loadings.

- col:

  `vector` of colours for score axes and loading axes and points/texts.

- expand:

  `numeric` expansion for the scores, defaulting to 1.

- labels:

  optional. If `"names"`, row names are used as labels. If `"numbers"`,
  row numbers are used as labels. (Can also be a vector of labels.)

- legendpos:

  `character` position of legend.

- ...:

  Additional arguments to `plot` and `scoreplot`.

## Value

No return, only a plot.

## Examples

``` r
# Load candies data
data(candies)

# Basic ASCA model with two factors and interaction
mod <- asca(assessment ~ candy * assessor, data=candies)

# Biplot
biplot(mod)


# Biplot with named loadings
biplot(mod, labels="names")

```
