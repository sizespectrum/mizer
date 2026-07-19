# Print mizer objects

Mizer supplies `print()` methods for the array-like objects returned by
many rate and summary functions. These methods print a compact preview
of the underlying matrix, array or vector: a header reporting the value
name, dimensions and units, followed by the actual values, truncated to
fit the console when the array is large. Species are truncated to a
leading subset, sizes to an evenly log-spaced sample spanning the full
size range (because size grids are uniform in log-space, this shows
small, medium and large individuals rather than just the smallest), and
time series to a representative sample of time steps that always
includes the first and last. A trailing note reports how much was
omitted. A three-dimensional
[`ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpeciesBySize.md)
object is previewed via its final time slice, matching the default
behaviour of
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) for that
class.

## Usage

``` r
# S3 method for class 'ArraySpeciesBySize'
print(x, ...)
# S3 method for class 'ArrayTimeBySpecies'
print(x, ...)
# S3 method for class 'ArrayTimeBySpeciesBySize'
print(x, ...)
# S3 method for class 'summary.ArraySpeciesBySize'
print(x, ...)
# S3 method for class 'summary.ArrayTimeBySpecies'
print(x, ...)
# S3 method for class 'summary.ArrayTimeBySpeciesBySize'
print(x, ...)
```

## Arguments

- x:

  The object to print.

- ...:

  Further arguments. They are currently ignored by the mizer methods.

## Value

The printed object, invisibly.

## Details

For full numeric access, use the object itself as an ordinary matrix,
array or vector, or convert it to a long data frame with
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md).

## See also

[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md),
[`ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md),
[`ArrayTimeBySpecies()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md),
[`ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpeciesBySize.md)

## Examples

``` r
# \donttest{
enc <- getEncounter(NS_params)
print(enc)
#> Encounter rate (12 species x 100 sizes) [g/year] 
#>          w
#> sp            0.001     0.083      5.78       480     39900
#>   Sprat   0.2992076  5.644995  86.66687  1415.973  39573.31
#>   Sandeel 0.4528175  8.595592 139.51701  2337.603  45507.81
#>   N.pout  0.5019776  8.287777 156.34311 13066.882 147886.62
#>   Herring 0.5752333 10.919341 177.48877  3117.693  80375.40
#>   Dab     0.4916095  8.431667 138.67700  7768.097 266704.99
#>   Whiting 0.4362525  7.245051 144.21893 14369.422 137539.48
#>   Sole    0.3646753  6.316326  95.11713  3313.102 148784.12
#>   Gurnard 0.3122260  5.367747  79.36674  3106.902 135351.68
#>   Plaice  0.2323659  3.947683  70.59929  4552.016 113242.16
#>   Haddock 0.5964130 10.429942 156.03134  5779.822 334719.58
#>   Cod     0.9658343 16.222960 277.59361 25639.777 436916.91
#>   Saithe  0.7709631 12.790324 174.21669  8311.746 187775.81
#> ... showing 5 of 100 sizes (0.001-39900 g, log-spaced); use as.data.frame() for the full data. 

biomass <- getBiomass(NS_sim)
print(biomass)
#> Biomass (44 times x 12 species) [g] 
#>       sp
#> time         Sprat      Sandeel       N.pout      Herring         Dab
#>   1967 50836384568 3.652370e+12 369838053933 524328920325 12000395984
#>   1973 50839649651 3.473508e+12 379087254642 275699594393 15382241615
#>   1979 53959238723 3.491434e+12 366843337845 387672220331 15243035268
#>   1985 26402174254 1.560356e+12 259806502396 470470157764 17635191965
#>   1992 23968401712 1.054672e+12 280632920267 430699336339 17284466852
#>   1998 31797035009 1.196296e+12 298030880150 429335887928 17479136749
#>   2004 39165801326 8.759191e+11 308055191214 466925281965 17886308035
#>   2010 36176429846 1.485299e+12 314880676670 408224168356 13458725479
#>       sp
#> time        Whiting         Sole      Gurnard       Plaice      Haddock
#>   1967 212534518393 128000031283 112908385702 2.276969e+12 627337734289
#>   1973 198345867124  98802757219 126999037442 2.504534e+12 540856360154
#>   1979 203662469253 103336430875 128644143259 2.405585e+12 549722419283
#>   1985 230237668108 101386030179 136973131775 2.250395e+12 660555395410
#>   1992 210367548378 119317446198 139092975786 2.113924e+12 635310276953
#>   1998 225880531337 100964781131 130801173560 1.974046e+12 686858953059
#>   2004 256540673031 103173835279 120019621650 2.104679e+12 797473701793
#>   2010 234404573351 104251576449 106732083088 2.756253e+12 820508475928
#>       sp
#> time            Cod       Saithe
#>   1967 1.618389e+12 1.062987e+12
#>   1973 4.168507e+11 8.263818e+11
#>   1979 3.953260e+11 7.936321e+11
#>   1985 3.145550e+11 6.917586e+11
#>   1992 3.253935e+11 7.026828e+11
#>   1998 2.851026e+11 8.633221e+11
#>   2004 2.624314e+11 9.193638e+11
#>   2010 4.036150e+11 7.748924e+11
#> ... showing 8 of 44 times (1967-2010, evenly spaced); use as.data.frame() for the full data. 
# }
```
