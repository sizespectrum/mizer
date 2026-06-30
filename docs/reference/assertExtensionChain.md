# Assert that an object's extension chain is compatible with the session

Stops with an informative error if the object's extension chain is not a
suffix of the session's registered maximal chain, or (when `check_class`
is `TRUE`) if the object does not inherit from the expected S4 marker
class.

## Usage

``` r
assertExtensionChain(
  object,
  extensions = objectExtensions(object),
  check_class = TRUE
)
```

## Arguments

- object:

  A `MizerParams` or `MizerSim` object.

- extensions:

  Named character vector giving the object's extension chain. Defaults
  to
  [`objectExtensions()`](https://sizespectrum.org/mizer/reference/objectExtensions.md)
  applied to `object`.

- check_class:

  Logical. If `TRUE` (default), also verify that `object` inherits from
  the expected S4 marker class.

## Value

Invisibly `TRUE`. Called for its side-effect of stopping on
incompatibility.
