# Get mean trophic level of each species

**\[experimental\]** Calculates the consumption-rate-weighted mean
trophic level of each species, defined as \$\$ T_i = \frac{\int
r_i(w)\\N_i(w)\\T_i(w)\\dw} {\int r_i(w)\\N_i(w)\\dw}, \$\$ where
\\r_i(w) = (1 - f_i(w))\\E_i(w)\\ is the consumption rate of an
individual of species \\i\\ at weight \\w\\, \\N_i(w)\\ is the abundance
density, and \\T_i(w)\\ is the size-resolved trophic level from
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md).
As in
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
the resource is given a size-dependent trophic level controlled by the
`w_R` and `beta_R` arguments.

## Usage

``` r
getTrophicLevelBySpecies(
  params,
  n = initialN(params),
  n_pp = initialNResource(params),
  n_other = initialNOther(params),
  w_R = 1e-10,
  beta_R = 1000,
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object.

- n:

  A matrix of species abundances (species x size). Defaults to the
  initial abundances stored in `params`.

- n_pp:

  A vector of the resource abundance by size. Defaults to the initial
  resource abundance stored in `params`.

- n_other:

  A named list of the abundances of other dynamical components. Defaults
  to the initial values stored in `params`.

- w_R:

  An average size (in grams) of primary producers in the resource
  spectrum, used to set the size-dependent resource trophic level.
  Defaults to `1e-10`.

- beta_R:

  An average predator/prey mass ratio for the resource spectrum, used to
  set the size-dependent resource trophic level. Must be greater than
  `1`. Defaults to `1000`.

- ...:

  Unused

## Value

A named vector with the mean trophic level for each species.

## See also

[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md)

Other summary functions:
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md),
[`getGrowthCurves()`](https://sizespectrum.org/mizer/reference/getGrowthCurves.md),
[`getN()`](https://sizespectrum.org/mizer/reference/getN.md),
[`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md),
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md),
[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)

## Examples

``` r
getTrophicLevelBySpecies(NS_params)
#>    Sprat  Sandeel   N.pout  Herring      Dab  Whiting     Sole  Gurnard 
#> 3.845388 3.702670 4.723848 3.596897 4.760884 4.958668 4.751892 4.302460 
#>   Plaice  Haddock      Cod   Saithe 
#> 4.417412 4.481840 5.241173 5.200990 
```
