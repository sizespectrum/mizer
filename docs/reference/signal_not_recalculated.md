# Signal that a rate array was not recalculated because it is frozen

Raised by the rate setters when they leave a frozen array alone although
the species parameters say that it should have a different value. It is
reported as a message, and `info_level = 0` silences it along with the
other information. Where no handler is collecting, for example when a
rate setter is called directly rather than via
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
it is shown anyway, because it may then be all the user hears. The
stronger
[`signal_frozen()`](https://sizespectrum.org/mizer/reference/signal_frozen.md)
warning is raised elsewhere, by whoever knows that the user asked for a
change, see
[`signal_frozen_changes()`](https://sizespectrum.org/mizer/reference/signal_frozen_changes.md).

## Usage

``` r
signal_not_recalculated(
  var,
  quantity,
  reset_call,
  derived_from = "species parameters"
)
```

## Arguments

- var:

  A string naming the slot that was not recalculated.

- quantity:

  A string naming the quantity for the user, for example "metabolic
  rate".

- reset_call:

  A string with the call that recalculates the quantity, for example
  "setMetabolicRate(params, reset = TRUE)".

- derived_from:

  A string naming the parameters that the quantity would have been
  calculated from.

## Value

`NULL` invisibly. Called for its side effect of signalling.

## Examples

``` r
with_info_level(
    signal_not_recalculated("metab", "metabolic rate",
                            "setMetabolicRate(params, reset = TRUE)")
)
#> The metabolic rate has been set manually and so is not recalculated from the species parameters. Call `setMetabolicRate(params, reset = TRUE)` if you want the metabolic rate to be calculated from the species parameters again.
```
