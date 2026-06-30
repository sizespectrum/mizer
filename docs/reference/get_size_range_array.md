# Get size range array

Helper function that returns an array (species x size) of logical values
indicating whether that size bin is within the size limits specified by
the arguments. Either the size limits can be the same for all species or
they can be specified as vectors with one value for each species in the
model.

## Usage

``` r
get_size_range_array(
  params,
  min_w = min(params@w),
  max_w = max(params@w),
  min_l = NULL,
  max_l = NULL,
  ...
)
```

## Arguments

- params:

  MizerParams object

- min_w:

  Smallest weight in size range. Defaults to smallest weight in the
  model.

- max_w:

  Largest weight in size range. Defaults to largest weight in the model.

- min_l:

  Smallest length in size range. If supplied, this takes precedence over
  `min_w`.

- max_l:

  Largest length in size range. If supplied, this takes precedence over
  `max_w`.

- ...:

  Unused

## Value

A logical array (species x size), with dimnames `sp` and `w`.

## Length to weight conversion

If `min_l` is specified there is no need to specify `min_w` and so on.
However, if a length is specified (minimum or maximum) then it is
necessary for the species parameter data.frame to include the parameters
`a` and `b` that determine the relation between length \\l\\ and weight
\\w\\ by \$\$w = a l^b.\$\$

It is possible to mix length and weight constraints, e.g. by supplying a
minimum weight and a maximum length, but this must be done the same for
all species. The default values are the minimum and maximum weights of
the spectrum, i.e., the full range of the size spectrum is used.
