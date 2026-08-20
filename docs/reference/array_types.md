# Kinds of quantity a mizer array can hold

Mizer arrays record what kind of quantity their values are in their
`type` attribute, because some kinds need handling that the numbers
alone do not reveal:

- `"value"`:

  the default: a rate, an amount, anything that needs no special
  handling.

- `"density"`:

  an amount per gram of body weight, like a number density. Plotting a
  density against a length axis restates it per centimetre, which
  changes the values and not just the axis.

- `"proportion"`:

  a fraction, like the feeding level. Plotted on a linear y axis showing
  the whole of the interval from 0 to 1, so that the value can be read
  against the scale it belongs to.

## Usage

``` r
array_types
```

## Format

A character vector of the three types.

## Details

A `"proportion"` is not *restricted* to the interval from 0 to 1: the
critical feeding level and the resource level can both exceed 1, and
their plots show it. The type is a statement about what the number
means, not a bound that mizer enforces.
