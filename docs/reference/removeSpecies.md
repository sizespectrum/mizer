# Remove species

This function removes all entries from the MizerParams object that refer
to the selected species. It does not recalculate the steady state for
the remaining species or retune their reproductive efficiency.

## Usage

``` r
removeSpecies(
  params,
  species,
  compensate = FALSE,
  info_level = default_info_level(),
  ...
)
```

## Arguments

- params:

  A mizer params object for the original system.

- species:

  The species to be removed. A vector of species names, or a numeric
  vector of species indices, or a logical vector indicating for each
  species whether it is to be removed (TRUE) or not.

- compensate:

  If `TRUE`, the external rates of the remaining species are adjusted
  and checked to preserve growth, diffusion, mortality and recruitment
  at the initial state. Supported resource dynamics are balanced. See
  the section "Compensating for the removed species". Default is
  `FALSE`.

- info_level:

  Controls the amount of information messages that are shown. Higher
  levels lead to more messages. Set to 0 to suppress all such messages.

- ...:

  Currently unused.

## Value

An object of type
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)

## Details

Removing species changes the rates of the species that remain: the
predators of the removed species lose some of their food, so their
encounter rate drops, and the prey of the removed species lose some of
their predators, so their predation mortality drops. A model that was at
steady state is therefore no longer at steady state after the removal.
If you want the remaining species not to notice the removal, set
`compensate = TRUE`, as described in the section "Compensating for the
removed species" below.

If a gear was targeting only the removed species, then this function
will NOT remove that gear. If you want to also remove that gear then you
can do that by calling
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md).

## Compensating for the removed species

With `compensate = TRUE`, the external rates of the remaining species
take over what the removed species contributed to their rates:

- the food that the remaining species found in the removed species is
  added to their external encounter rate, see
  [`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),

- the predation mortality that the removed species imposed on the
  remaining species is added to their external mortality rate, see
  [`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md),

- the diffusion that the removed prey species caused is added to the
  external diffusion rate, see
  [`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md).
  This only matters if
  [`use_predation_diffusion()`](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)
  is switched on.

Compensation is evaluated through the complete effective rate pipeline
at the initial state. After trimming species and running all declared
hooks, mizer checks growth, diffusion, total mortality and final
recruitment of every remaining species. A failed check stops the
operation. These quantities are preserved at this state; their responses
to later changes may differ.

For the built-in semichemostat and logistic resource, reduced predation
is balanced analytically by changing capacity, using the effective
resource mortality. This puts the resource at steady state and keeps its
renewal rate. A capacity set manually is replaced. A constant resource
needs no adjustment. Custom resource dynamics and other components
require explicit removal support; their hooks receive
`details$compensate` and are responsible for their own compensation.
This is not a whole-model steady-state guarantee.

The core adds each external-rate shortfall once. Custom callbacks that
scale or otherwise alter that response need an authored transformation
hook with permission to update the affected arrays. Mizer does not
estimate a numerical correction for an unknown callback.

Each external rate that is changed is marked as set manually. A later
change to the species parameters that the rate is usually calculated
from (such as `E_ext`, `z0` or `D_ext`) therefore does not overwrite the
compensation; mizer tells you when it ignores such a change. External
rates to which the removed species contributed nothing are left
unchanged and unmarked.

## See also

[`addSpecies()`](https://sizespectrum.org/mizer/reference/addSpecies.md),
[`renameSpecies()`](https://sizespectrum.org/mizer/reference/renameSpecies.md)

## Examples

``` r
params <- NS_params
species_params(params)$species
#>  [1] "Sprat"   "Sandeel" "N.pout"  "Herring" "Dab"     "Whiting" "Sole"   
#>  [8] "Gurnard" "Plaice"  "Haddock" "Cod"     "Saithe" 
params <- removeSpecies(params, c("Cod", "Haddock"))
species_params(params)$species
#>  [1] "Sprat"   "Sandeel" "N.pout"  "Herring" "Dab"     "Whiting" "Sole"   
#>  [8] "Gurnard" "Plaice"  "Saithe" 

# Removing Cod and Haddock moves the model off its steady state, unless the
# remaining species are compensated for their loss
rowSums(getSteadyResidual(removeSpecies(NS_params, c("Cod", "Haddock"))))
#>       Sprat     Sandeel      N.pout     Herring         Dab     Whiting 
#>  0.24068904  0.09865627  0.45749567  0.42999736  0.38410262  0.42358061 
#>        Sole     Gurnard      Plaice      Saithe 
#>  0.25579948  0.43635238  0.30499751 -0.01867903 
rowSums(getSteadyResidual(removeSpecies(NS_params, c("Cod", "Haddock"),
                                        compensate = TRUE)))
#>       Sprat     Sandeel      N.pout     Herring         Dab     Whiting 
#> 0.003278363 0.001320613 0.013855625 0.008077633 0.005050589 0.003729652 
#>        Sole     Gurnard      Plaice      Saithe 
#> 0.003286390 0.009183216 0.006807205 0.002111153 
```
