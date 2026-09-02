# Reconcile the species parameters with the given species parameters

**\[experimental\]** Records among the given species parameters every
value in `species_params(params)` that mizer would otherwise change when
it recalculates the species parameters, so that the model's species
parameters become a fixed point: no recalculation, however often
repeated, moves them again. Called automatically by
[`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md).

## Usage

``` r
reconcileSpeciesParams(params, info_level = default_info_level())
```

## Arguments

- params:

  A MizerParams object.

- info_level:

  Controls the amount of information messages and warnings that are
  shown. Higher levels lead to more messages, `info_level = 0` gives
  silence. The default is taken from the `mizer_info_level` option, see
  [`default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.md).

## Value

The MizerParams object with the additional given species parameters
recorded.

## Details

Mizer distinguishes between the species parameters that were given
explicitly and the ones it calculated itself, see
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md).
Only the given ones are protected: whenever the species parameters are
recalculated – which every use of `species_params<-()` triggers – the
calculated ones are derived afresh from the given ones.

A value that was written straight into the `species_params` slot, with
for example `params@species_params$h[1] <- 20`, is therefore not
protected and the next parameter change undoes it without saying so.
Models saved before mizer kept track of the given species parameters,
and models built by code that changes the slot directly without also
calling
[`record_given_species_params()`](https://sizespectrum.org/mizer/reference/record_given_species_params.md),
can hold many such values.

This function repairs such a model. It works out what the species
parameters would look like if they were recalculated from the given
species parameters now, compares that against the species parameters the
model actually holds, and records the entries that differ among the
given species parameters. It then repeats that with the enlarged record,
and keeps going until a recalculation reproduces the species parameters
exactly.

The repetition is what makes the result a fixed point. Recording one
value generally changes what a recalculation gives for the parameters
mizer derives from it, so a single pass is not enough. If you had set
`h` by hand, for instance, the model's `gamma` was derived from the old
`h`; the first pass records `h`, and only the second notices that
`gamma` would now be derived afresh from your new `h` and records it
too. Once the loop finishes, nothing in
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)
moves again, no matter how many times the species parameters are
recalculated.

The comparison is made entry by entry, so a parameter is recorded only
for the species whose value differs, and only where the model actually
holds a value: an `NA` in
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)
is not a value the user put there and is never recorded. A value that a
recalculation already reproduces is left calculated and goes on
responding to changes in the parameters it is derived from.

The model is not changed by this function: no species parameter value,
and no rate array, is touched. Only the record of where the values came
from is.

## The price of the fixed point

Freezing the model as it stands means recording values that mizer
calculated rather than values you supplied – the `gamma` of the example
is mizer's own, derived from an `h` that is no longer in the model. That
is deliberate: it is the `gamma` the model's rate arrays were built
from, and the alternative is to let a load or an unrelated parameter
change silently alter the model.

The cost is that such a parameter no longer responds to the parameters
it was derived from. To hand one back to mizer's calculation, clear its
entry to `NA` in `given_species_params(params)`, as for any other given
species parameter.

## See also

[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md),
[`given_species_params()`](https://sizespectrum.org/mizer/reference/species_params.md),
[`record_given_species_params()`](https://sizespectrum.org/mizer/reference/record_given_species_params.md)

## Examples

``` r
params <- NS_params
# Write a value straight into the species parameter slot, bypassing the
# record of the values the user has supplied.
params@species_params$w_mat25[1] <- species_params(params)$w_mat25[1] / 2
# Mizer does not know that this value was given, so a recalculation would
# undo it.
is.null(given_species_params(params)$w_mat25)
#> [1] TRUE

params <- reconcileSpeciesParams(params)
#> The species parameter `w_mat25` holds a value that a recalculation would not reproduce. I have recorded it among the given species parameters so that it is not overwritten.
# Now it is recorded, for the one species whose value differs.
given_species_params(params)$w_mat25
#>   Sprat Sandeel  N.pout Herring     Dab Whiting    Sole Gurnard  Plaice Haddock 
#> 5.82373      NA      NA      NA      NA      NA      NA      NA      NA      NA 
#>     Cod  Saithe 
#>      NA      NA 
```
