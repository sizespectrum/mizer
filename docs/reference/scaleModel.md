# Change scale of the model

**\[experimental\]**

The abundances in mizer and some rates depend on the size of the area to
which they refer. So they could be given per square meter or per square
kilometre or for an entire study area or any other choice of yours. This
function allows you to change the scale of the model by automatically
changing the abundances and rates accordingly.

## Usage

``` r
scaleModel(params, factor, ...)
```

## Arguments

- params:

  A MizerParams object

- factor:

  The factor by which the scale is multiplied

- ...:

  Additional arguments passed to the method.

## Value

The rescaled MizerParams object

## Details

If you rescale the model by a factor \\c\\ then this function makes the
following rescalings in the params object:

- The initial consumer and resource abundances are rescaled by \\c\\.

- The search volume is rescaled by \\1/c\\.

- The resource carrying capacity is rescaled by \\c\\

- Reproduction rate parameters `R_max` and `constant_reproduction` are
  multiplied by \\c\\; `ricker_b` and `sheperd_b` are divided by it.

The effect of this is that the dynamics of the rescaled model are
identical to those of the unscaled model, in the sense that it does not
matter whether one first calls `scaleModel()` and then runs a simulation
with [`project()`](https://sizespectrum.org/mizer/reference/project.md)
or whether one first runs a simulation and then rescales the resulting
abundances.

Components and custom providers must declare transformation support.
Their hooks run after the core fields are scaled and before final
validation, with the original model available for reference. A component
can scale a stock while leaving intensive quantities, such as
temperature, unchanged. Mizer does not infer how to transform an
arbitrary component state.

In practice you will need to use some observations to set the scale for
your model. If you have biomass observations you can use
[`calibrateBiomass()`](https://sizespectrum.org/mizer/reference/calibrateBiomass.md),
if you have observed numbers you can use
[`calibrateNumber()`](https://sizespectrum.org/mizer/reference/calibrateNumber.md).
