# Summarise mizer objects

Mizer provides `summary()` methods for model objects and for the
specialised array classes returned by many mizer functions.

## Usage

``` r
# S3 method for class 'ArraySpeciesBySize'
summary(object, all.sizes = FALSE, ...)
# S3 method for class 'ArrayTimeBySpecies'
summary(object, ...)
# S3 method for class 'ArrayTimeBySpeciesBySize'
summary(object, all.sizes = FALSE, ...)
# S3 method for class 'MizerSim'
summary(object, ...)
# S3 method for class 'MizerParams'
summary(object, ...)
```

## Arguments

- object:

  The object to summarise.

- all.sizes:

  If `FALSE` (the default), values outside a species' size range
  (`w_min` to `w_max`) are left out, as in
  [`plot()`](https://sizespectrum.org/mizer/reference/plot.md). Only for
  the classes with a size dimension.

- ...:

  Further arguments. They are currently ignored by the mizer methods.

## Value

For
[`MizerParams()`](https://sizespectrum.org/mizer/reference/MizerParams.md)
and
[`MizerSim()`](https://sizespectrum.org/mizer/reference/MizerSim.md),
the object is returned invisibly. For array objects, a list of class
`summary.ArraySpeciesBySize`, `summary.ArrayTimeBySpecies` or
`summary.ArrayTimeBySpeciesBySize`.

## Details

For a
[`MizerParams()`](https://sizespectrum.org/mizer/reference/MizerParams.md)
object, `summary()` prints the model metadata, size grids, selected
species parameters and fishing gear details. For a
[`MizerSim()`](https://sizespectrum.org/mizer/reference/MizerSim.md)
object, it first prints the parameter summary and then reports the
simulated time period and output interval.

For
[`ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md),
[`ArrayTimeBySpecies()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md)
and
[`ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpeciesBySize.md)
objects, `summary()` returns a small list with the value name, units,
dimensions and a per-species data frame containing minimum, mean and
maximum values. Printing that summary object gives the same compact
table in a human-readable form.

For the two classes that have a size dimension, those values are taken
over each species' own size range, from its `w_min` to its `w_max`,
which is the range
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) draws. A
rate array is defined on the whole size grid, but the values outside a
species' range describe an animal that does not exist — the encounter
rate a 40 kg Sprat would have — and they are usually the extreme ones,
so a summary that included them reported the size grid rather than the
species. Pass `all.sizes = TRUE` for the whole grid.

## See also

[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md),
[`str()`](https://sizespectrum.org/mizer/reference/str.md),
[`MizerParams()`](https://sizespectrum.org/mizer/reference/MizerParams.md),
[`MizerSim()`](https://sizespectrum.org/mizer/reference/MizerSim.md),
[`ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md),
[`ArrayTimeBySpecies()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md),
[`ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpeciesBySize.md)

## Examples

``` r
# \donttest{
summary(NS_params)
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> An object of class "MizerParams" 
#> mizer version: 3.0.0.9003
#> Created: 2021-09-03 20:29:38
#> Modified: 2026-06-24 13:50:03
#> Consumer size spectrum:
#>  minimum size:   0.001
#>  maximum size:   39851.3
#>  no. size bins:  100
#> Resource size spectrum:
#>  minimum size:   2.12182e-13
#>  maximum size:   9.82091
#>  no. size bins:  179 (226 size bins in total)
#> Steady state:
#>  biomass drift:  0.014 /year (at steady state)
#> Species details:
#> An object of class "species_params" containing parameters for 12 species:
#>  species   w_inf w_mat w_min  f0   beta sigma
#>    Sprat    33.0    13 0.001 0.6  51076   0.8
#>  Sandeel    36.0     4 0.001 0.6 398849   1.9
#>   N.pout   100.0    23 0.001 0.6     22   1.5
#>  Herring   334.0    99 0.001 0.6 280540   3.2
#>      Dab   324.0    21 0.001 0.6    191   1.9
#>  Whiting  1192.0    75 0.001 0.6     22   1.5
#>     Sole   866.0    78 0.001 0.6    381   1.9
#>  Gurnard   668.0    39 0.001 0.6    283   1.8
#>   Plaice  2976.0   105 0.001 0.6    113   1.6
#>  Haddock  4316.5   165 0.001 0.6    558   2.1
#>      Cod 39851.3  1606 0.001 0.6     66   1.3
#>   Saithe 39658.6  1076 0.001 0.6     40   1.1
#> 
#> Fishing gear details:
#> Gear          Effort  Target species 
#>  ----------------------------------
#> Industrial     0.00   Sprat, Sandeel, N.pout 
#> Pelagic        1.00   Herring 
#> Beam           0.50   Dab, Sole, Plaice 
#> Otter          0.50   Whiting, Gurnard, Haddock, Cod, Saithe 
summary(NS_sim)
#> An object of class "MizerSim" 
#> Parameters:
#> An object of class "MizerParams" 
#> mizer version: 3.0.0.9003
#> Created: 2021-09-03 20:30:02
#> Modified: 2026-06-24 13:50:03
#> Consumer size spectrum:
#>  minimum size:   0.001
#>  maximum size:   39851.3
#>  no. size bins:  100
#> Resource size spectrum:
#>  minimum size:   8.72744e-13
#>  maximum size:   9.82091
#>  no. size bins:  171 (218 size bins in total)
#> Steady state:
#>  biomass drift:  0.91 /year  (not at steady state, largest in Cod - run tuneSteadyState())
#> Species details:
#> An object of class "species_params" containing parameters for 12 species:
#>  species   w_inf w_mat w_min   beta sigma
#>    Sprat    33.0    13 0.001  51076   0.8
#>  Sandeel    36.0     4 0.001 398849   1.9
#>   N.pout   100.0    23 0.001     22   1.5
#>  Herring   334.0    99 0.001 280540   3.2
#>      Dab   324.0    21 0.001    191   1.9
#>  Whiting  1192.0    75 0.001     22   1.5
#>     Sole   866.0    78 0.001    381   1.9
#>  Gurnard   668.0    39 0.001    283   1.8
#>   Plaice  2976.0   105 0.001    113   1.6
#>  Haddock  4316.5   165 0.001    558   2.1
#>      Cod 39851.3  1606 0.001     66   1.3
#>   Saithe 39658.6  1076 0.001     40   1.1
#> 
#> Fishing gear details:
#> Gear          Effort  Target species 
#>  ----------------------------------
#> Sprat          0.51   Sprat 
#> Sandeel        0.56   Sandeel 
#> N.pout         0.51   N.pout 
#> Herring        1.29   Herring 
#> Dab            0.93   Dab 
#> Whiting        0.73   Whiting 
#> Sole           1.14   Sole 
#> Gurnard        0.28   Gurnard 
#> Plaice         0.93   Plaice 
#> Haddock        0.68   Haddock 
#> Cod            0.94   Cod 
#> Saithe         0.74   Saithe 
#>  Note: effort varied over time for Sprat (0.00 to 1.76), Sandeel (0.00 to 1.48), N.pout (0.00 to 1.76), Herring (0.12 to 3.33), Dab (0.40 to 1.40), Whiting (0.30 to 1.07), Sole (0.76 to 1.58), Gurnard (0.00 to 1.00), Plaice (0.40 to 1.40), Haddock (0.18 to 1.08), Cod (0.68 to 1.10), Saithe (0.33 to 1.36); mean shown above.
#> Simulation parameters:
#>  Time period: 1967 to 2010
#>  Output stored every 1 years
#>  Time step   Method: 
summary(getEncounter(NS_params))
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> Encounter rate [g/year] 
#> 12 species x 100 sizes
#> 
#>  Species       Min        Mean         Max
#>    Sprat 0.2992076    37.92979    239.6706
#>  Sandeel 0.4528175    67.08030    434.6202
#>   N.pout 0.5019776   237.84164   2613.4140
#>  Herring 0.5752333   282.49597   2203.3821
#>      Dab 0.4916095   395.38383   4473.2748
#>  Whiting 0.4362525  2507.48128  28974.8976
#>     Sole 0.3646753   474.92455   5616.2042
#>  Gurnard 0.3122260   318.45335   3726.8753
#>   Plaice 0.2323659  1898.07284  24800.9828
#>  Haddock 0.5964130  3547.12149  49555.6495
#>      Cod 0.9658343 52646.61006 436916.9135
#>   Saithe 0.7709631 14646.02285 160506.0739
summary(getFMort(NS_sim))
#> Fishing mortality [1/year] 
#> 44 times x 12 species x 100 sizes
#> 
#>  Species          Min       Mean       Max
#>    Sprat 0.0000000000 0.12879552 2.1827923
#>  Sandeel 0.0000000000 0.11611087 1.3076336
#>   N.pout 0.0000000000 0.13718398 2.1228994
#>  Herring 0.0077894362 0.18281500 1.4419293
#>      Dab 0.0017746971 0.02322755 0.1633433
#>  Whiting 0.0126246154 0.17579846 1.3138868
#>     Sole 0.0268320133 0.15939358 1.0259979
#>  Gurnard 0.0000000000 0.00541050 0.1054407
#>   Plaice 0.0089964121 0.18786555 0.8670560
#>  Haddock 0.0015854377 0.22672527 1.4275500
#>      Cod 0.0494762790 0.36512704 1.0721072
#>   Saithe 0.0009633361 0.15208058 1.2032857
# }
```
