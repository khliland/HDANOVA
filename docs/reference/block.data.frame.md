# Block-wise indexable data.frame

This is a convenience function for making `data.frame`s that are easily
indexed on a block-wise basis.

## Usage

``` r
block.data.frame(X, block_inds = NULL, to.matrix = TRUE)
```

## Arguments

- X:

  Either a single `data.frame` to index or a `list` of
  matrices/data.frames

- block_inds:

  Named `list` of indexes if `X` is a single `data.frame`, otherwise
  `NULL`.

- to.matrix:

  `logical` indicating if input list elements should be converted to
  matrices.

## Value

A `data.frame` which can be indexed block-wise.

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
[`asca_plots`](https://khliland.github.io/HDANOVA/reference/asca_plots.md),
[`pcanova_results`](https://khliland.github.io/HDANOVA/reference/pcanova_results.md)
and
[`pcanova_plots`](https://khliland.github.io/HDANOVA/reference/pcanova_plots.md)

## Examples

``` r
# Random data
M <- matrix(rnorm(200), nrow = 10)
# .. with dimnames
dimnames(M) <- list(LETTERS[1:10], as.character(1:20))

# A named list for indexing
inds <- list(B1 = 1:10, B2 = 11:20)

X <- block.data.frame(M, inds)
str(X)
#> 'data.frame':    10 obs. of  2 variables:
#>  $ B1: 'AsIs' num [1:10, 1:10] -1.8273 -0.3382 1.0092 -1.7455 -0.0482 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:10] "1" "2" "3" "4" ...
#>  $ B2: 'AsIs' num [1:10, 1:10] -0.821 2.083 0.53 -0.841 -1.249 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:10] "11" "12" "13" "14" ...
```
