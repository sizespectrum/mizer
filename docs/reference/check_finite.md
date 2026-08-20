# Check that the rate arrays hold only finite values

Gives an error if any of the arrays in the MizerParams object contains
non-finite values, with the exception of the maximum intake rate, which
is allowed to be infinite.

## Usage

``` r
check_finite(params)
```

## Arguments

- params:

  A MizerParams object.

## Value

`TRUE`, invisibly.

## Details

This check is cheap and is run on every call to
[`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md),
also when the repair work is skipped, because it catches exactly the
kind of damage that the fingerprint in
[`validation_key()`](https://sizespectrum.org/mizer/reference/validation_key.md)
does not cover: a bad value written into an array whose shape is
unchanged.
