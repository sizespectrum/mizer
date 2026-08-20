# Density measures a spectrum can be expressed in

A size spectrum is a density, and a density only has a meaning together
with the variable it is a density with respect to. That variable is one
of

- `"w"`:

  a density with respect to weight, e.g. numbers per gram.

- `"log_w"`:

  a density with respect to logarithmic weight, e.g. numbers per log
  weight interval.

- `"l"`:

  a density with respect to length, e.g. numbers per cm.

- `"log_l"`:

  a density with respect to logarithmic length.

- `NA`:

  not a density, e.g. a rate or a dimensionless quantity. Such values
  are left alone when the size axis changes.

## Usage

``` r
density_measures
```

## Format

A character vector of the four density measures.

## Details

Mizer arrays are indexed by the model's weight grid, so an array that
holds a density (`type = "density"`, see
[array_types](https://sizespectrum.org/mizer/reference/array_types.md))
always holds one with respect to weight. The other measures arise for
quantities the spectrum plots compute on the fly:
`plotSpectra(per_log_size = TRUE)` shows a density with respect to
logarithmic weight, and either can be restated per unit length by
`size_axis = "l"`.
