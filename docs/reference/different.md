# Check whether two objects are different

Check whether two objects are numerically different, ignoring all
attributes.

## Usage

``` r
different(a, b)
```

## Arguments

- a:

  First object

- b:

  Second object

## Value

TRUE or FALSE

## Details

We use this helper function in particular to see if a new value for a
slot in MizerParams is different from the existing value in order to
give the appropriate messages.
