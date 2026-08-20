# Expand kernel weights indexed by grid offset into a full kernel array

The predation kernel depends only on the predator/prey mass ratio, so on
the geometric grid it is a function of the offset \\m\\ between the
predator and prey grid indices alone. `phis[i, m + 1]` holds the weight
of species \\i\\ at offset \\m\\. This helper writes those weights into
the (predator species x predator size x prey size) array that the
non-FFT code paths work with.

## Usage

``` r
expand_kernel_offsets(phis, params, species)
```

## Arguments

- phis:

  A species-by-offset matrix of kernel weights, with the offset running
  from 0 to `length(params@w_full) - 1`.

- params:

  A MizerParams object supplying the grid.

- species:

  A character vector of species names for the dimnames.

## Value

An array (predator species x predator size x prey size).
