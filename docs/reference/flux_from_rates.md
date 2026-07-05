# Assemble the flux matrix from growth, diffusion and recruitment rates

Internal helper holding the arithmetic shared by `getFlux.MizerParams`
and `getFlux.MizerSim`. Keeping it separate lets the `MizerSim` method
resolve the rate functions once and reuse them across all saved time
steps.

## Usage

``` r
flux_from_rates(params, n, g, d, rdd, power = 0, flux_limiter = "none")
```

## Arguments

- params:

  A valid `MizerParams` object.

- n:

  A matrix of species abundances (species x size).

- g:

  Growth rate matrix (species x size), as from
  [`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md).

- d:

  Diffusion rate matrix (species x size), as from
  [`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md).

- rdd:

  Density-dependent reproduction rate vector (one per species), as from
  [`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md).

- power:

  The flux at weight \\w\\ is multiplied by \\w\\ raised to `power`. The
  default `power = 0` leaves the flux of individuals unchanged.

- flux_limiter:

  Advective-flux scheme: `"none"` (first-order upwind), `"van_leer"` or
  `"centred"` (second-order log-size scheme). Defaults to `"none"`.

## Value

A plain species x size matrix of fluxes (no mizer array class).
