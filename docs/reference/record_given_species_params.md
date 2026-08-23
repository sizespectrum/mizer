# Record the species parameters that have changed

**\[experimental\]** Compares the new species parameters in `value`
against the old ones in `old_sp` and records the entries that have
changed in the given species parameter data frame `given`. This is the
change detection used by `species_params<-()`, exported so that code
which updates the species parameters by other means can record its
changes the same way.

## Usage

``` r
record_given_species_params(given, value, old_sp)
```

## Arguments

- given:

  The given species parameter data frame to record into, usually
  `given_species_params(params)`.

- value:

  A data frame holding the new species parameters.

- old_sp:

  A data frame holding the species parameters as they were before the
  change. Must have one row per species, in the same order as `value`
  and `given`.

## Value

The updated `given` data frame.

## Details

Mizer distinguishes between the species parameters that were given
explicitly and those that it calculated itself, see
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md).
Only the given ones are protected: whenever the species parameters are
recalculated – which every use of `species_params<-()` triggers – the
calculated ones are derived afresh from the given ones. So code that
computes a species parameter and writes it into the `species_params`
slot directly has its work silently undone by the next parameter change,
unless the value is also recorded among the given species parameters.

The usual way to record a parameter is to set it with
`species_params<-()`, which also rebuilds the species parameters and
recalculates all the rates that depend on them. This function is the
recording step on its own, for the case where the caller has already
updated the affected rates itself, for example an optimiser that fits a
species parameter and the rate array it determines together. Rebuilding
and recalculating would then be wasted work, and can even undo the
caller's own adjustment.

Only the values that have actually changed are recorded. This matters:
recording an unchanged value would turn a calculated species parameter
into a given one and thereby stop it from responding to changes in the
parameters it is derived from. The comparison is made entry by entry, so
a parameter is protected only for the species whose value changed. `NA`
is compared as a value rather than as an unknown, so `NA` staying `NA`
does not count as a change. A column that is not present in `old_sp` at
all is taken to be new and is recorded in full.

## See also

[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md),
[`given_species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)

## Examples

``` r
params <- NS_params
sp_before <- species_params(params)
given_before <- given_species_params(params)

# Set a species parameter and the rate it determines, without going through
# `species_params<-()` and its recalculation of every other rate.
params@species_params$ks[1] <- species_params(params)$ks[1] * 2
params@metab[1, ] <- params@metab[1, ] * 2

# Record the change so that it is not recalculated away later
params@given_species_params <-
    record_given_species_params(given_species_params(params),
                                species_params(params), sp_before)

# Only the entry that changed has been recorded
given_species_params(params)$ks == given_before$ks
#>   Sprat Sandeel  N.pout Herring     Dab Whiting    Sole Gurnard  Plaice Haddock 
#>   FALSE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
#>     Cod  Saithe 
#>    TRUE    TRUE 
```
