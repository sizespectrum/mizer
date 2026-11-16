# Keep resource abundance constant

If you set your resource dynamics to use this function then the resource
abundances are kept constant over time.

## Usage

``` r
resource_constant(params, n_pp, ...)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.md)
  object

- n_pp:

  A vector of the resource abundance by size

- ...:

  Unused

## Value

Vector containing the resource number density in each size class at the
next timestep

## Details

To set your model to keep the resource constant over time you do

    resource_dynamics(params) <- "resource_constant"

where you should replace `params` with the name of the variable holding
your MizerParams object.

## See also

[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md)

Other resource dynamics functions:
[`resource_logistic()`](https://sizespectrum.org/mizer/reference/resource_logistic.md),
[`resource_semichemostat()`](https://sizespectrum.org/mizer/reference/resource_semichemostat.md)

## Examples

``` r
params <- NS_params
resource_dynamics(params) <- "resource_constant"
```
