# ASCA Plot Methods

Various plotting procedures for
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.md) objects.

## Usage

``` r
# S3 method for class 'asca'
loadingplot(object, factor = 1, comps = 1:2, ...)

# S3 method for class 'asca'
scoreplot(
  object,
  factor = 1,
  comps = 1:2,
  within_level = "all",
  pch.scores = 19,
  pch.projections = 1,
  gr.col = NULL,
  projections = TRUE,
  spider = FALSE,
  ellipsoids,
  confidence,
  xlim,
  ylim,
  xlab,
  ylab,
  legendpos,
  ...
)

permutationplot(object, factor = 1, xlim, xlab = "SSQ", main, ...)

rotationplot(object, factor = 1, xlim, xlab = "SSQ", main, ...)
```

## Arguments

- object:

  `asca` object.

- factor:

  `integer/character` for selecting a model factor. If factor \<= 0 or
  "global", the PCA of the input is used (negativ factor to include
  factor level colouring with global PCA).

- comps:

  `integer` vector of selected components.

- ...:

  additional arguments to underlying methods.

- within_level:

  MSCA parameter for chosing plot level (default = "all").

- pch.scores:

  `integer` plotting symbol.

- pch.projections:

  `integer` plotting symbol.

- gr.col:

  `integer` vector of colours for groups.

- projections:

  Include backprojections in score plot (default = TRUE).

- spider:

  Draw lines between group centers and backprojections (default =
  FALSE).

- ellipsoids:

  `character` "confidence" or "data" ellipsoids for balanced fixed
  effect models.

- confidence:

  `numeric` vector of ellipsoid confidences, default = c(0.4, 0.68,
  0.95).

- xlim:

  `numeric` x limits.

- ylim:

  `numeric` y limits.

- xlab:

  `character` x label.

- ylab:

  `character` y label.

- legendpos:

  `character` position of legend.

- main:

  Plot title.

## Value

The plotting routines have no return.

## Details

Usage of the functions are shown using generics in the examples in
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.md). Plot
routines are available as `scoreplot.asca` and `loadingplot.asca`.

## References

- Smilde, A., Jansen, J., Hoefsloot, H., Lamers,R., Van Der Greef, J.,
  and Timmerman, M.(2005). ANOVA-Simultaneous Component Analysis (ASCA):
  A new tool for analyzing designed metabolomics data. Bioinformatics,
  21(13), 3043–3048.

- Liland, K.H., Smilde, A., Marini, F., and Næs,T. (2018). Confidence
  ellipsoids for ASCA models based on multivariate regression theory.
  Journal of Chemometrics, 32(e2990), 1–13.

- Martin, M. and Govaerts, B. (2020). LiMM-PCA: Combining ASCA+ and
  linear mixed models to analyse high-dimensional designed data. Journal
  of Chemometrics, 34(6), e3232.

## See also

Main methods:
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.md),
[`apca`](https://khliland.github.io/HDANOVA/reference/apca.md),
[`limmpca`](https://khliland.github.io/HDANOVA/reference/limmpca.md),
[`msca`](https://khliland.github.io/HDANOVA/reference/msca.md),
[`pcanova`](https://khliland.github.io/HDANOVA/reference/pcanova.md),
[`prc`](https://khliland.github.io/HDANOVA/reference/prc.md) and
[`permanova`](https://khliland.github.io/HDANOVA/reference/permanova.md).
Workhorse function underpinning most methods:
[`hdanova`](https://khliland.github.io/HDANOVA/reference/hdanova.md).
Extraction of results and plotting:
[`asca_results`](https://khliland.github.io/HDANOVA/reference/asca_results.md),
`asca_plots`,
[`pcanova_results`](https://khliland.github.io/HDANOVA/reference/pcanova_results.md)
and
[`pcanova_plots`](https://khliland.github.io/HDANOVA/reference/pcanova_plots.md)
