# Principal Response Curves

Wrapper for the
[`prc`](https://vegandevs.github.io/vegan/reference/prc.html) function
to allow for formula input.

## Usage

``` r
prc(formula, data, ...)
```

## Arguments

- formula:

  Model formula accepting a single response (block) and predictors. If
  no predictor is called 'time', time is assumed to be the second
  predictor.

- data:

  The data set to analyse.

- ...:

  Additional arguments to
  [`prc`](https://vegandevs.github.io/vegan/reference/prc.html).

## Value

An object of class `prc`.

## See also

Main methods:
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.md),
[`apca`](https://khliland.github.io/HDANOVA/reference/apca.md),
[`limmpca`](https://khliland.github.io/HDANOVA/reference/limmpca.md),
[`msca`](https://khliland.github.io/HDANOVA/reference/msca.md),
[`pcanova`](https://khliland.github.io/HDANOVA/reference/pcanova.md),
`prc` and
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
data(caldana)
(pr <- prc(compounds ~ light * time, caldana))
#> 
#> Call: prc(formula = compounds ~ light * time, data = caldana)
#> 
#>                Inertia Proportion Rank
#> Total         11.47745    1.00000     
#> Conditional    1.11209    0.09689    6
#> Constrained    2.51542    0.21916   18
#> Unconstrained  7.84993    0.68394   67
#> 
#> Inertia is variance
#> 
#> Eigenvalues for constrained axes:
#>   RDA1   RDA2   RDA3   RDA4   RDA5   RDA6   RDA7   RDA8   RDA9  RDA10  RDA11 
#> 0.9779 0.4049 0.3239 0.1755 0.1540 0.1247 0.0699 0.0607 0.0498 0.0336 0.0311 
#>  RDA12  RDA13  RDA14  RDA15  RDA16  RDA17  RDA18 
#> 0.0273 0.0245 0.0187 0.0168 0.0094 0.0066 0.0061 
#> 
#> Eigenvalues for unconstrained axes:
#>    PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
#> 1.8824 1.4836 0.7107 0.4964 0.3533 0.3435 0.3159 0.2450 
#> (Showing 8 of 67 unconstrained eigenvalues)
#> 
summary(pr)
#> 
#> Call:
#> prc(formula = compounds ~ light * time, data = caldana) 
#> Species scores:
#>                    Alanine                     Valine 
#>                  -0.900828                  -0.175768 
#>                    Leucine                 Isoleucine 
#>                   0.657559                   0.264687 
#>                    Proline                     Serine 
#>                  -1.141926                   0.051640 
#>                  Threonine               beta-alanine 
#>                  -0.075174                   0.726471 
#>             Hydroxyproline                       GABA 
#>                   0.007045                   1.087892 
#>                  Aspartate                 Asparagine 
#>                   0.696076                   0.338664 
#>                 Methionine            O-acetyl-serine 
#>                  -0.863337                  -0.205739 
#>                  Glutamate              Phenylalanine 
#>                   0.337782                  -1.245032 
#>                  Ornithine                  Glutamine 
#>                  -0.041141                  -1.302687 
#>                     Lysine                   Tyrosine 
#>                   1.175389                   0.606966 
#>              Threonic-acid        Citrulline-Arginine 
#>                   0.323093                  -0.526770 
#>               Pyruvic-acid                Citric-acid 
#>                  -0.473488                   0.414141 
#>              Succinic-acid               Fumaric-acid 
#>                   0.690860                  -0.005060 
#>                 Malic-acid                Lactic-acid 
#>                  -0.006521                   0.497413 
#>              Glycolic-acid               Benzoic-acid 
#>                  -1.624649                   0.141741 
#>                Maleic-acid             Nicotinic-acid 
#>                  -0.767982                   0.060987 
#>              Itaconic-acid                Citramalate 
#>                   0.169012                  -0.050705 
#>     4-hydroxy-benzoic-acid Dehydroascorbic-acid-dimer 
#>                   0.191369                  -0.062572 
#>              Gluconic-acid       Dehydroascorbic-acid 
#>                   0.159411                   0.108530 
#>              Ascorbic-acid     4-Hydroxycinnamic-acid 
#>                   0.232603                   0.133527 
#>         Similar-to-Adenine                  Shikimate 
#>                  -0.010357                  -0.878436 
#>                 Erythritol                  Arabinose 
#>                   0.068112                   0.110114 
#>                   Arabitol                     Fucose 
#>                   0.337455                   0.051972 
#>                   Fructose                   Mannitol 
#>                  -3.606930                  -0.215745 
#>                  Galactose                    Glucose 
#>                  -0.097331                  -2.777645 
#>                    Sucrose                    Maltose 
#>                  -0.615215                   0.263945 
#>                  Trehalose                 Galactinol 
#>                   0.195746                  -0.211351 
#>               myo-inositol                     Uracil 
#>                   0.070161                   1.067653 
#>                 Putrescine               Ethanolamine 
#>                   0.210320                   0.546516 
#>                   Glycerol      Indole-3-acetonitrile 
#>                   0.146868                  -0.185760 
#>               Sinapic-acid              Palmitic-acid 
#>                  -0.131682                   0.267378 
#>          Octadecanoic-acid            Docosanoic-acid 
#>                   0.361142                   0.299480 
#>         Tetracosanoic-acid          Hexacosanoic-acid 
#>                  -0.060571                  -0.188997 
#>          Octacosanoic-acid 
#>                  -0.167153 
#> 
#> Coefficients for treatment + time:treatment interaction
#> which are contrasts to treatment Dark 
#> rows are treatment, columns are time
#>                     0        5       10        20        40      80     160
#> Low Light   2.465e-16 0.003685  0.05458  0.121789  0.006526 -0.1098 -0.1833
#> Light      -3.826e-17 0.038250 -0.10345  0.084547 -0.098440 -0.2416 -0.3368
#> High Light  2.047e-16 0.010944  0.01016 -0.009301 -0.312173 -0.6746 -0.7716
```
