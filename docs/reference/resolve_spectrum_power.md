# Resolve the power of weight multiplying a spectrum

The quantity plotted by
[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md)
is the number density multiplied by `w^power`. That power is the sum of
two independent choices: whether the quantity is a biomass density (a
factor of `w`) or a number density, and whether it is a density with
respect to logarithmic size (another factor of `w`) or with respect to
size. The two choices are what determine the y-axis label and the
Jacobian used when converting to a length axis, and they are not
recoverable from `power` alone: `power = 1` is both the biomass density
with respect to weight and the number density with respect to
logarithmic weight.

## Usage

``` r
resolve_spectrum_power(power = NULL, biomass = NULL, per_log_size = NULL)
```

## Arguments

- power:

  The power of weight multiplying the number density, or `NULL`.

- biomass:

  Whether to plot a biomass density rather than a number density, or
  `NULL`. The default is `TRUE`.

- per_log_size:

  Whether to plot a density with respect to logarithmic size rather than
  with respect to size, or `NULL`. The default is `FALSE`.

## Value

A list with entries `power`, `biomass` and `per_log_size`.

## Details

Users therefore express the choice with the `biomass` and `per_log_size`
flags. The `power` argument remains available, both for backwards
compatibility and as an escape hatch for powers that are not the sum of
two flags. Each argument is `NULL` when it was not supplied by the user.
