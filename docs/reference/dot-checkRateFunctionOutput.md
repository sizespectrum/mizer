# Check that a rate function returns the correct output dimensions

Called by
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
to verify that a candidate rate function returns an array (or
vector/list) of the correct dimensions for the requested rate.

## Usage

``` r
.checkRateFunctionOutput(params, rate, fun)
```

## Arguments

- params:

  A MizerParams object

- rate:

  Name of the rate being replaced, e.g. `"Encounter"`.

- fun:

  Name of the candidate function to validate.

## Value

Invisibly `NULL`. Called for its side-effect of stopping with an
informative error if the output has the wrong shape.
