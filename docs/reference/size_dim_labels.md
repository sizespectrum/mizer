# Identify the dimensions of an array over the size grid

Internal helper for
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md).
Returns a label for each dimension of `x`. The last dimension must run
over the size grid and is labelled `"w"`. Other dimensions are labelled
from the names of their dimnames, if they have any, with `"species"` and
`"size"` normalised to mizer's `"sp"` and `"w"`. An unnamed
second-to-last dimension is labelled `"sp"` if its extent is the number
of species. Any remaining unnamed dimension gets a unique label of its
own, so that it is carried through to the result rather than matched
against a dimension of the abundance.

## Usage

``` r
size_dim_labels(x, arg, no_sp, no_w)
```

## Arguments

- x:

  The array to label.

- arg:

  The name of the argument holding `x`, used in error messages.

- no_sp:

  The number of species in the model.

- no_w:

  The number of size bins in the model.

## Value

A character vector with one label for each dimension of `x`.

## Details

A scalar has no dimensions and gets no labels.
