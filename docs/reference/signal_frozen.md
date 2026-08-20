# Signal that a change the user made cannot take effect

A rate array that has been set by hand is protected by a comment, see
the "Setting or changing rates" section in
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md).
Mizer then no longer calculates it from the species parameters, so a
change to one of the species parameters that feeds it has no effect on
the model. This function raises the condition that tells the user so.

## Usage

``` r
signal_frozen(var, message)
```

## Arguments

- var:

  A string naming the quantity the report is about.

- message:

  The message to give the user.

## Value

`NULL` invisibly. Called for its side effect of signalling.

## Details

The condition is raised at severity `"warning"`, see
[`signal_info()`](https://sizespectrum.org/mizer/reference/signal_info.md),
so that it survives the
[`suppressMessages()`](https://rdrr.io/r/base/message.html) that
[`species_params<-()`](https://sizespectrum.org/mizer/reference/species_params.md)
runs over its recalculation. It also carries the class
`info_about_frozen` for code that wants to catch this kind of report in
particular.

Only signal this when the user has actually asked for something that is
not happening. The mere fact that a frozen array differs from what the
formula would give is not enough: mizer freezes arrays itself when it
builds the trait-based and community models, and those arrays differ
from the formula for the lifetime of the model. See
[`signal_frozen_changes()`](https://sizespectrum.org/mizer/reference/signal_frozen_changes.md),
which decides this from the species parameters the user changed.
