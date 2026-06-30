# Get or set the use_predation_diffusion flag

Controls whether predation-induced diffusion is included when
calculating rates with
[`mizerDiffusion()`](https://sizespectrum.org/mizer/reference/mizerDiffusion.md).
When `FALSE` (the default), the predation-driven diffusion term is
omitted, preserving the behaviour of previous mizer versions. Set to
`TRUE` to enable the diffusion term from the jump-growth equation.

## Usage

``` r
use_predation_diffusion(params)

use_predation_diffusion(params) <- value
```

## Arguments

- params:

  A MizerParams object.

- value:

  A single logical value (`TRUE` or `FALSE`).

## Value

`use_predation_diffusion()`: A single logical value.

`use_predation_diffusion<-`: A MizerParams object with the
`use_predation_diffusion` flag updated.

## See also

Other functions for setting parameters:
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md),
[`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md),
[`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),
[`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md),
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md),
[`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.md),
[`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md),
[`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md),
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md),
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md),
[`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md),
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)
