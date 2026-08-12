# Mizer

> An R package for dynamic multi-species size-spectrum modelling of fish communities.

Mizer models marine ecosystems subject to fishing. It tracks the full size distribution of each species and the plankton resource, computing growth, predation, and mortality rates from individual-level physiology and feeding interactions. A model is set up from a species parameter data frame, calibrated to a steady state matching observations, then projected forward under different fishing strategies.

## Installation

```r
install.packages("mizer")
remotes::install_github("sizespectrum/mizer")  # development version
```

## Quick start

```r
library(mizer)
params <- newMultispeciesParams(NS_species_params, NS_interaction)
sim <- project(params, t_max = 10, effort = 0)
plot(sim)
```

See the [Get started](https://sizespectrum.org/mizer/articles/mizer.md) article for a full walkthrough. Extension packages add new biology (temperature dependence, starvation mortality, seasonal dynamics) — see [Using extension packages](https://sizespectrum.org/mizer/articles/using-extension-packages.md).
