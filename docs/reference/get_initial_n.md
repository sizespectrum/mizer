# Calculate initial population abundances

This function uses the model parameters and other parameters to
calculate initial values for the species number densities. These initial
abundances are currently quite arbitrary and not close to the steady
state. We intend to improve this in the future.

## Usage

``` r
get_initial_n(params, n0_mult = NULL, a = 0.35)
```

## Arguments

- params:

  The model parameters. An object of type
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md).

- n0_mult:

  Multiplier for the abundance at size 0 when using defaults edition 1.
  If not supplied, `kappa / 1000` is used. This argument is ignored for
  defaults edition 2 and later.

- a:

  A parameter with a default value of 0.35.

## Value

An `ArraySpeciesBySize` object (species x size) of population
abundances.

## Examples

``` r
init_n <- get_initial_n(NS_params)
```
