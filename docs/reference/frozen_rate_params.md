# Which parameters feed which frozen array

A lookup table used by
[`signal_frozen_changes()`](https://sizespectrum.org/mizer/reference/signal_frozen_changes.md)
to decide whether a change the user made can take effect. Each entry is
named after a slot of
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
that can be frozen and gives the quantity as the user knows it, the call
that unfreezes it, the parameters that the setter and its default
calculations read, and what kind of parameters those are.

## Usage

``` r
frozen_rate_params()
```

## Value

A named list of lists with entries `quantity`, `reset_call`, `params`
and `derived_from`.

## Details

The list of parameters does not have to be exhaustive, and deliberately
is not: it names the parameters that the setter reads directly together
with the main inputs of the default calculations for those parameters. A
parameter that is missing simply means that the user is not warned, and
is left with the message that the setter itself gives, see
[`signal_not_recalculated()`](https://sizespectrum.org/mizer/reference/signal_not_recalculated.md).
Listing a parameter that in fact has no influence is the worse mistake,
because it warns about a change that did take effect.
