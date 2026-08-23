# The default level of information that mizer gives

**\[experimental\]** Returns the `mizer_info_level` option if it is set
and `fallback` otherwise. This is the default of the `info_level`
argument of the functions that report information, so that
`options(mizer_info_level = 0)` quietens mizer as a whole, including the
functions that have no `info_level` argument of their own, such as
[`species_params<-()`](https://sizespectrum.org/mizer/reference/species_params.md)
and the rate setters.

## Usage

``` r
default_info_level(fallback = 3)
```

## Arguments

- fallback:

  The level to use when the option is not set. Defaults to 3, which
  reports everything.

## Value

A single number, or `NA` to leave the reporting to a handler further
out.

## Details

Extension packages should use this as the default of their own
`info_level` argument, so that a constructor or setter of theirs follows
the option like mizer's own do:

    newFooParams <- function(species_params, ...,
                             info_level = default_info_level()) {
        newMultispeciesParams(species_params, info_level = info_level, ...)
    }

Take the argument explicitly like that rather than hard-coding a value
in the call, which would make a user's own `info_level` collide with it.

## Examples

``` r
default_info_level()
#> [1] 3

# Setting the option changes what every reporting function defaults to.
old <- options(mizer_info_level = 1)
default_info_level()
#> [1] 1
options(old)
```
