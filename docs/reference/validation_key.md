# Fingerprint of the slots that determine the outcome of the repair and structural validity checks in `validParams()`

The fingerprint covers every slot that
[`repair_params()`](https://sizespectrum.org/mizer/reference/repair_params.md)
reads or writes and every slot that the `MizerParams` validity function
inspects, except for the values inside the large rate arrays, of which
only the dimensions and dimension names are included. The values in
those arrays are checked unconditionally by
[`check_finite()`](https://sizespectrum.org/mizer/reference/check_finite.md)
instead.

## Usage

``` r
validation_key(params)
```

## Arguments

- params:

  A MizerParams object.

## Value

A string.

## Details

The fingerprint is always calculated afresh from the current contents of
the object, never stored on the object itself, so no change to the
object can escape it.
