# Test whether one extension chain is a suffix of another

An empty `candidate` is always a suffix. Order and values must match
exactly for the overlapping tail.

## Usage

``` r
isSuffixChain(candidate, chain)
```

## Arguments

- candidate:

  Named character vector to test.

- chain:

  Named character vector that may contain `candidate` as a tail.

## Value

`TRUE` if `candidate` is a suffix of `chain`, `FALSE` otherwise.
