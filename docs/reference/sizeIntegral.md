# Integrate a quantity over the size spectrum

**\[experimental\]** Calculates \$\$\int\_{w\_{min}}^{w\_{max}} N_i(w)\\
K_i(w)\\ dw\$\$ for each species \\i\\, using the quadrature scheme that
the model is actually using. This is the recommended way to write your
own summary or indicator function: it selects the size range, applies
the bin-averaging appropriate to the model's
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)
setting and wraps the result in the appropriate mizer array class, so
that none of those rules need to be remembered. The built-in summary
functions like
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
[`getN()`](https://sizespectrum.org/mizer/reference/getN.md),
[`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md) and
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md) are
all implemented with it.

## Usage

``` r
sizeIntegral(
  object,
  weighting = 1,
  n = NULL,
  ...,
  value_name = NULL,
  units = NULL
)
```

## Arguments

- object:

  A `MizerParams` or a `MizerSim` object.

- weighting:

  The weighting factor \\K(w)\\ of the integral, evaluated on the size
  grid. See the section "The weighting factor" below. Defaults to 1.

- n:

  The abundance density. Either a species x size matrix or a time x
  species x size array. Defaults to the initial abundance
  `initialN(object)` for a `MizerParams` object and to the saved
  abundances `object@n` for a `MizerSim` object.

- ...:

  Arguments passed to
  [`get_size_range_array()`](https://sizespectrum.org/mizer/reference/get_size_range_array.md)
  to select the size range to integrate over, i.e. `min_w`, `max_w`,
  `min_l` and `max_l`.

- value_name:

  A string giving a human-readable name for the value, used when the
  result is wrapped in a mizer array class.

- units:

  A string giving the units of the result, used when the result is
  wrapped in a mizer array class.

## Value

The value of the integral, see the section "Shape of the result" above.

## The weighting factor

The weighting factor \\K(w)\\ is supplied already evaluated on the size
grid. It can be

- a single number (the default `weighting = 1` integrates the abundance
  density itself, giving numbers),

- a vector with one value for each size bin, which is then used for all
  species,

- a matrix (species x size), for example `params@maturity`,

- an array with further dimensions in front, for example the gear x
  species x size array returned by
  [`getFMortGear()`](https://sizespectrum.org/mizer/reference/getFMortGear.md)
  or the time x species x size array returned by `getFMort(sim)`. Those
  extra dimensions are carried through to the result.

If the weighting factor is a product of several size-dependent factors,
pass the whole product: bin-averaging is applied to the product as a
single weighting factor, which is not the same as averaging the factors
separately.

Do **not** include the bin widths `params@dw` in the weighting factor
and do not bin-average it yourself; `sizeIntegral()` does both.

## Shape of the result

The size dimension is integrated out. The remaining dimensions are those
of `n` together with any extra dimensions of `weighting`, so

- with a `MizerParams` object and a species x size weighting the result
  is a named vector with one value per species,

- with a `MizerSim` object it is an
  [ArrayTimeBySpecies](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md)
  object (time x species),

- with a gear x species x size weighting the extra `gear` dimension is
  kept, giving a gear x species array (or time x gear x species for a
  `MizerSim`).

Dimensions of `weighting` other than the last two are matched to the
dimensions of `n` by the names of their dimnames, so a weighting whose
first dimension is named `"time"` is lined up with the times of the
simulation rather than producing an outer product.

## See also

[`get_size_range_array()`](https://sizespectrum.org/mizer/reference/get_size_range_array.md),
[`bin_average_weight()`](https://sizespectrum.org/mizer/reference/bin_average_weight.md),
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)

Other summary functions:
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md),
[`getGrowthCurves()`](https://sizespectrum.org/mizer/reference/getGrowthCurves.md),
[`getN()`](https://sizespectrum.org/mizer/reference/getN.md),
[`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md),
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md),
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
[`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md),
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md),
[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)

## Examples

``` r
# The biomass of each species, i.e. what getBiomass() does
sizeIntegral(NS_params, weighting = NS_params@w)
#>        Sprat      Sandeel       N.pout      Herring          Dab      Whiting 
#> 4.054293e+11 5.589441e+12 4.762654e+11 1.467576e+12 1.056360e+10 1.965198e+11 
#>         Sole      Gurnard       Plaice      Haddock          Cod       Saithe 
#> 1.433825e+11 6.237753e+10 1.404752e+12 4.581486e+11 6.015918e+11 5.213786e+11 

# ... restricted to a size range
sizeIntegral(NS_params, weighting = NS_params@w, min_w = 10, max_w = 1000)
#>        Sprat      Sandeel       N.pout      Herring          Dab      Whiting 
#> 2.441192e+11 4.589606e+12 2.388476e+11 1.273446e+12 8.373096e+09 1.615158e+11 
#>         Sole      Gurnard       Plaice      Haddock          Cod       Saithe 
#> 1.274005e+11 2.475989e+10 7.660875e+11 3.331725e+11 4.515965e+10 1.589695e+11 

# The numbers of individuals larger than 10g
sizeIntegral(NS_params, min_w = 10)
#>        Sprat      Sandeel       N.pout      Herring          Dab      Whiting 
#>  13565641567 203862451898   7823764384  26562731929    240601035   2308726264 
#>         Sole      Gurnard       Plaice      Haddock          Cod       Saithe 
#>   2512241983    962788847  16641583382   5953056353    412647877   2604150904 

# Spawning stock biomass: the weighting is the product maturity * w
K <- sweep(NS_params@maturity, 2, NS_params@w, "*")
sizeIntegral(NS_params, weighting = K)
#>        Sprat      Sandeel       N.pout      Herring          Dab      Whiting 
#> 2.108102e+11 5.378411e+12 1.831597e+11 4.426440e+11 6.885676e+09 1.135728e+11 
#>         Sole      Gurnard       Plaice      Haddock          Cod       Saithe 
#> 6.363024e+10 9.102233e+09 3.033659e+11 1.519652e+11 5.315367e+11 3.275782e+11 

# An indicator through time, ready to plot
biomass <- sizeIntegral(NS_sim, weighting = NS_params@w,
                        value_name = "Biomass", units = "g")
biomass[c("1972", "2010"), c("Herring", "Cod")]
#> Biomass (2 times x 2 species) [g] 
#>       sp
#> time        Herring          Cod
#>   1972 218218354800 436573537805
#>   2010 408224168356 403614998737
```
