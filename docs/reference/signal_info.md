# Signal information about a choice mizer made

**\[experimental\]** Raises the condition that
[`with_info_level()`](https://sizespectrum.org/mizer/reference/with_info_level.md)
collects. This is the way for mizer, and for anything extending it, to
tell the user about a default it filled in, an input it adjusted or an
instruction it could not carry out, without deciding on its own how
loudly to say it: the handler installed by whichever function the user
actually called does that.

## Usage

``` r
signal_info(
  var,
  message,
  level = 3,
  severity = c("info", "warning"),
  unhandled = c("drop", "show"),
  class = character()
)
```

## Arguments

- var:

  A string naming the quantity the report is about.

- message:

  The message to give the user.

- level:

  How important the report is. Level 1 is important enough to survive
  `info_level = 1`, level 3 is chatter that only the default
  `info_level = 3` shows.

- severity:

  `"info"` to report as a message, `"warning"` to report as a warning.
  Use `"warning"` when the user asked for something that is not
  happening, because a message can be, and on the
  [`species_params<-()`](https://sizespectrum.org/mizer/reference/species_params.md)
  path is, suppressed.

- unhandled:

  What to do when no handler is collecting, for example because a rate
  setter was called directly rather than through
  [`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md).
  `"drop"` says nothing, which suits chatter that only makes sense as
  part of a report about a whole model. `"show"` reports it there and
  then, at the same severity: a message for `"info"` and a warning for
  `"warning"`.

- class:

  Further classes to give the condition, for code that wants to catch a
  particular kind of report.

## Value

`NULL` invisibly. Called for its side effect of signalling.

## Details

Progress reports are the one thing that does not belong here: they have
to appear while the work is going on, and these are collected and given
at the end.

## Examples

``` r
# With nothing collecting, a `"drop"` report says nothing at all ...
signal_info("h", "Using a default for `h`.")

# ... whereas `unhandled = "show"` reports it there and then.
signal_info("h", "Using a default for `h`.", unhandled = "show")
#> Using a default for `h`.

# Normally it is raised inside a call whose body is wrapped in
# `with_info_level()`, which is what decides whether to show it.
with_info_level(signal_info("h", "Using a default for `h`."))
#> Using a default for `h`.
```
