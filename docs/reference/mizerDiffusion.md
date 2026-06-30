# Calculate diffusion rate

Calculates the diffusion rate \\D_i(w)\\ (grams^2/year) for each
species. This diffusion rate has two components:

1.  The diffusion due due to the variability in prey sizes. This is the
    diffusion term from the jump-growth equation.

2.  Any externally specified diffusion, which is added via
    [`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md)

You would not usually call this function directly but instead use
[`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md),
which then calls this function unless an alternative diffusion rate
function has been registered, see
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md).

## Usage

``` r
projectDiffusion(params, n, n_pp, n_other, t = 0, feeding_level, ...)

# S3 method for class 'MizerParams'
projectDiffusion(params, n, n_pp, n_other, t = 0, feeding_level, ...)

mizerDiffusion(params, n, n_pp, n_other, t = 0, feeding_level, ...)
```

## Arguments

- params:

  A MizerParams object

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size

- n_other:

  A list of abundances for other dynamical components

- t:

  The time for which to do the calculation (Not used by standard mizer
  rate functions but useful for extensions.)

- feeding_level:

  An array (species x size) with the feeding level. If not provided, it
  is calculated from the given abundances.

- ...:

  Unused

## Value

A two dimensional array (species x size) holding the diffusion rate.
