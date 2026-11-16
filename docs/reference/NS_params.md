# Example MizerParams object for the North Sea example

A MizerParams object created from the `NS_species_params_gears` species
parameters and the `inter` interaction matrix together with an initial
condition corresponding to the steady state obtained from fishing with
an effort
`effort = c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 0.5)`.

## Usage

``` r
NS_params
```

## Format

A MizerParams object

## Source

Blanchard et al.

## See also

Other example parameter objects:
[`NS_sim`](https://sizespectrum.org/mizer/reference/NS_sim.md)

## Examples

``` r
# \donttest{
sim = project(NS_params, effort = c(Industrial = 0, Pelagic = 1, 
                                    Beam = 0.5, Otter = 0.5))
plot(sim)

# }
```
