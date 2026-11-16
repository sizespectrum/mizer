# Initial values for fish spectra

Values used as starting values for simulations with
[`project()`](https://sizespectrum.org/mizer/reference/project.md).

## Usage

``` r
initialN(params) <- value

initialN(object)
```

## Arguments

- params:

  A MizerParams object

- value:

  A matrix with dimensions species x size holding the initial number
  densities for the fish spectra.

- object:

  An object of class MizerParams or MizerSim

## Value

An `ArraySpeciesBySize` object with dimensions species x size holding
the initial number densities for the fish spectra.

## See also

[`initialNResource()`](https://sizespectrum.org/mizer/reference/initialNResource-set.md),
[`initialNOther()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md)

## Examples

``` r
# Doubling abundance of Cod in the initial state of the North Sea model
params <- NS_params
initialN(params)["Cod", ] <- 2 * initialN(params)["Cod", ]
```
