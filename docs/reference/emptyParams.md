# Create empty MizerParams object of the right size

An internal function. Sets up a valid
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
object with all the slots initialised and given dimension names, but
with some slots left empty. This function is to be used by other
functions to set up full parameter objects.

## Usage

``` r
emptyParams(
  species_params,
  gear_params = data.frame(),
  no_w = 100,
  min_w = 0.001,
  max_w = NA,
  min_w_pp = 1e-12
)
```

## Arguments

- species_params:

  A data frame of species-specific parameter values.

- gear_params:

  A data frame with gear-specific parameter values.

- no_w:

  The number of size bins in the consumer spectrum.

- min_w:

  Sets the size of the eggs of all species for which this is not given
  in the `w_min` column of the `species_params` dataframe.

- max_w:

  The largest size of the consumer spectrum. By default this is set to
  the largest `w_max` specified in the `species_params` data frame.

- min_w_pp:

  The smallest size of the resource spectrum.

## Value

An empty but valid MizerParams object

## Size grid

A size grid is created so that the log-sizes are equally spaced. The
spacing is chosen so that there will be `no_w` fish size bins, with the
smallest starting at `min_w` and the largest starting at `max_w`. For
the resource spectrum there is a larger set of bins containing
additional bins below `min_w`, with the same log size. The number of
extra bins is such that `min_w_pp` comes to lie within the smallest bin.

## Changes to species params

The `species_params` slot of the returned MizerParams object may differ
from the data frame supplied as argument to this function because
default values are set for missing parameters.

## See also

See
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
for a function that fills the slots left empty by this function.
