# Initial value for resource spectrum

Value used as starting value for simulations with
[`project()`](https://sizespectrum.org/mizer/reference/project.md).

## Usage

``` r
initialNResource(params) <- value

initialNResource(object)
```

## Arguments

- params:

  A MizerParams object

- value:

  A vector with the initial number densities for the resource spectrum

- object:

  An object of class MizerParams or MizerSim

## Value

A vector with the initial number densities for the resource spectrum

## See also

[`initialN()`](https://sizespectrum.org/mizer/reference/initialN-set.md),
[`initialNOther()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md)

## Examples

``` r
# Doubling resource abundance in the initial state of the North Sea model
params <- NS_params
initialNResource(params) <- 2 * initialNResource(params)
```
