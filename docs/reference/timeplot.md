# Timeplot for Combined Effects

Timeplot for Combined Effects

## Usage

``` r
timeplot(
  object,
  factor,
  time,
  comb,
  comp = 1,
  ylim,
  x_time = FALSE,
  xlab = time,
  ylab = paste0("Score ", comp),
  lwd = 2,
  ...
)
```

## Arguments

- object:

  `asca` object.

- factor:

  `integer/character` main factor.

- time:

  `integer/character` time factor.

- comb:

  `integer/character` combined effect factor.

- comp:

  `integer` component number.

- ylim:

  `numeric` y limits.

- x_time:

  `logical` use time levels as non-equispaced x axis (default = FALSE).

- xlab:

  `character` x label.

- ylab:

  `character` y label.

- lwd:

  `numeric` line width.

- ...:

  additional arguments to `plot`.

## Value

Nothing

## Examples

``` r
data("caldana")
mod.comb <- asca(compounds ~ time + comb(light + time:light), data=caldana)

# Default time axis
timeplot(mod.comb, factor="light", time="time", comb=2)


# Non-equispaced time axis (using time levels)
timeplot(mod.comb, factor="light", time="time", comb=2, x_time=TRUE)


# Second component
timeplot(mod.comb, factor="light", time="time", comb=2, comp=2, x_time=TRUE)
```
