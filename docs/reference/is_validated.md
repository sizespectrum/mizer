# Keep track of which MizerParams objects have been fully validated

Because the fingerprint returned by
[`validation_key()`](https://sizespectrum.org/mizer/reference/validation_key.md)
determines the outcome of
[`repair_params()`](https://sizespectrum.org/mizer/reference/repair_params.md)
and of the structural validity checks, an object whose fingerprint has
been recorded by a previous call to
[`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md)
needs neither. The record lasts for the R session only.

## Usage

``` r
is_validated(key)

record_validated(key, max_size = validated_params_max)

clear_validated_params()
```

## Arguments

- key:

  A fingerprint as returned by
  [`validation_key()`](https://sizespectrum.org/mizer/reference/validation_key.md).

- max_size:

  The number of fingerprints to keep. When the record has grown to this
  size it is emptied.

## Value

`is_validated()` returns TRUE if the fingerprint has been recorded.

`record_validated()` returns `NULL`, invisibly.

`clear_validated_params()` returns `NULL`, invisibly. The next
validation of any object then takes the full path again.
