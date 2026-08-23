# Pick out some of the series of a scan

The series of a scan need not be species: a scan of
[`getMeanWeight()`](https://sizespectrum.org/mizer/reference/getMeanWeight.md)
has a single series that no model has ever heard of. So the selection is
made against the series the scan actually holds rather than through
[`valid_species_arg()`](https://sizespectrum.org/mizer/reference/valid_species_arg.md),
which would reject anything that is not a species in the model.

## Usage

``` r
select_scan_series(available, species)
```

## Arguments

- available:

  The `Species` column of the scan.

- species:

  The series asked for, as names, as whole-number indices into the
  series of the scan, or as a logical vector with one entry per series.
  Must select at least one.

## Value

A logical vector selecting the rows to keep.
