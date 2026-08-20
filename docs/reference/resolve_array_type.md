# Resolve the type of a mizer array

Called by the array constructors. An explicit `type` is validated and
used as given; `NULL` means the constructor was called without the
argument, in which case a density is recognised from the other metadata,
the way mizer recognised one before the `type` attribute existed. That
keeps arrays built by extension packages, and arrays saved by earlier
versions, behaving as they did.

## Usage

``` r
resolve_array_type(type, value_name = NULL, units = NULL)
```

## Arguments

- type:

  The type supplied to the constructor, or `NULL`.

- value_name:

  The `value_name` of the array.

- units:

  The `units` of the array.

## Value

One of
[array_types](https://sizespectrum.org/mizer/reference/array_types.md).
