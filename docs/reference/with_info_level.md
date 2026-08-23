# Collect and report the information signals raised while setting parameters

**\[experimental\]** While mizer sets up or changes a model it raises
conditions of class `info_about_default` to tell the user about the
choices it made on their behalf and about the instructions it could not
carry out. This function evaluates `expr` with a calling handler that
collects those conditions and reports them together once `expr` has
finished, so that the user gets one report rather than a stream of
messages.

## Usage

``` r
with_info_level(expr, info_level = default_info_level(), except = character())
```

## Arguments

- expr:

  The expression to evaluate. It is evaluated in the calling
  environment, so assignments made in it have the same effect as they
  would have without this wrapper.

- info_level:

  The level of information to report, or `NA` to leave the reporting to
  a handler further out. Defaults to
  [`default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.md),
  which consults the `mizer_info_level` option.

- except:

  A character vector of `var`s not to report on, for a caller that is
  going to say the same thing itself. Everything else raised inside
  `expr` is reported as usual.

## Value

The value of `expr`.

## Details

Each condition carries three fields, see
[`signal_info()`](https://sizespectrum.org/mizer/reference/signal_info.md):

- `var` names the quantity the report is about.

- `level` says how important it is, a low level meaning important and a
  high level meaning chatter. Only conditions with `level` at most
  `info_level` are reported.

- `severity` says how to report it: `"info"` conditions become a single
  [`message()`](https://rdrr.io/r/base/message.html) and `"warning"`
  conditions a single
  [`warning()`](https://rdrr.io/r/base/warning.html).

The severity matters because
[`species_params<-()`](https://sizespectrum.org/mizer/reference/species_params.md)
runs [`suppressMessages()`](https://rdrr.io/r/base/message.html) over
its recalculation to quieten the routine chatter. A report that the user
needs to see even there — that an instruction of theirs had no effect —
must therefore be a warning, see
[`signal_frozen()`](https://sizespectrum.org/mizer/reference/signal_frozen.md).

Identical reports are collapsed, so a quantity that is reported on twice
in the same call takes up one line, but two different things said about
the same quantity are both kept.

## Nesting

Handlers nest by themselves: while one is collecting, any handler
installed further in steps aside and lets the outer one do the
reporting. A function can therefore wrap its body in `with_info_level()`
without knowing whether its caller has already done so, which is what
allows every entry point to install a handler. `info_level = NA` asks
for the same thing explicitly, for the rare case where a function wants
to leave the reporting to a caller that has not installed a handler yet.

The reporting happens on exit, so a function can wrap its whole body
even though it returns from the middle of it.

Silence is the exception to "the outermost handler decides":
`info_level = 0` drops the reports raised inside it even when a handler
further out is collecting, so that a function can build something
quietly as part of a larger job that does report.

## Examples

``` r
# Wrap the body of a function that reports, and everything raised inside it
# is collected and given together once the call has finished.
myConstructor <- function(x, info_level = default_info_level()) {
    with_info_level(info_level = info_level, {
        signal_info("h", "No `h` provided, using a default.", level = 1)
        signal_info("gamma", "Calculating `gamma` from `f0`.")
        x
    })
}
myConstructor(1)
#> ℹ No `h` provided, using a default.
#> ℹ Calculating `gamma` from `f0`.
#> [1] 1

# `info_level = 1` keeps only the report that was marked important.
myConstructor(1, info_level = 1)
#> No `h` provided, using a default.
#> [1] 1

# `info_level = 0` is silence.
myConstructor(1, info_level = 0)
#> [1] 1
```
