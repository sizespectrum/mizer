# Restate the units of a density in a different density measure

The per-size factor is found in either of the two spellings mizer uses
for it, `1/g` and `g^-1`, and is then swapped for the size unit of the
target measure. A density per *logarithmic* size carries no size unit at
all — a number per log weight interval is just a number — so converting
to one removes the factor instead: `1/g` becomes dimensionless and
`g^-1/year` becomes `1/year`.

## Usage

``` r
convert_density_units(units, from, to)
```

## Arguments

- units:

  The units of the values, possibly `NULL`.

- from, to:

  Density measures, see
  [density_measures](https://sizespectrum.org/mizer/reference/density_measures.md).

## Value

The units expressed in the `to` measure.

## Details

Units that state no per-size factor are returned unchanged, since there
is then nothing to identify as the size unit.
