# Repair a MizerParams object

Rebuilds the parameter tables and the slots that are derived from them.
This is a deterministic function of a small number of slots and is
idempotent: on an already-repaired object it recomputes identical
values.
[`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md)
therefore skips it for an object whose fingerprint has already been
recorded, see
[`validation_key()`](https://sizespectrum.org/mizer/reference/validation_key.md)
and
[`is_validated()`](https://sizespectrum.org/mizer/reference/is_validated.md).

## Usage

``` r
repair_params(params)
```

## Arguments

- params:

  A MizerParams object.

## Value

The repaired MizerParams object.
