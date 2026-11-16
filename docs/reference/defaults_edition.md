# Default editions

Function to set and get which edition of default choices is being used.

## Usage

``` r
defaults_edition(edition = NULL)

defaults_edition(params) <- value
```

## Arguments

- edition:

  NULL, a number, or a MizerParams or MizerSim object.

- params:

  A MizerParams object.

- value:

  The edition to record in the model.

## Value

If `edition` is `NULL`, the currently active edition number. If
`edition` is a MizerParams or MizerSim object, the edition recorded in
that model. If `edition` is a number, the function sets the global
`mizer_defaults_edition` option, emits a message, and returns the
supplied value invisibly.

`defaults_edition<-`: The MizerParams object with the new edition
recorded. Nothing is recalculated.

## Details

The mizer functions for creating new models make a lot of choices for
default values for parameters that are not provided by the user.
Sometimes we find better ways to choose the defaults and update mizer
accordingly. When we do this, we will increase the edition number.

If you call `defaults_edition()` without an argument it returns the
currently active edition. If you call it with a number it sets the
active edition to that value.

Users who want their existing code for creating models not to change
behaviour when run with future versions of mizer should explicitly set
the desired defaults edition at the top of their code.

## Each model keeps its own edition

The active edition is used when a new model is created, and the model
records it. Whenever mizer later fills in a default for that model, for
example when it recalculates the model after you have changed a species
parameter, it uses the edition recorded in the model, not the one that
is active at the time. So changing the active edition affects only
models created afterwards, and an existing model keeps behaving as it
did when it was created. Call `defaults_edition(params)` to see the
edition of a model. Models created with a version of mizer that did not
yet record the edition are of edition 1.

To give an existing model the defaults of a different edition, create it
afresh with that edition active.

`defaults_edition(params) <- value` changes the edition recorded in a
model. It does not recalculate anything: the model keeps the values it
has, and the new edition is used only for defaults that mizer fills in
from then on. It is meant for a model that was created with a different
edition active, before mizer recorded the edition in the model, and has
therefore been marked as edition 1.

## The editions

The most recent edition is edition 3. The current default is edition 1.
Edition 2 will become the default in the next release. The following
defaults are changed in edition 2:

- `catchability` = 0.3 instead of 1

- `initial effort` = 1 instead of 0

- `gamma` is set to ensure a feeding level of `f0` for larvae with the
  current value of `interaction_resource` instead of
  interaction_resource = 1\`.

- `initial_n` is set using
  [`get_steady_state_n()`](https://sizespectrum.org/mizer/reference/get_steady_state_n.md)
  instead of the rather arbitrary old choice.

- In
  [`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md),
  `psi` is no longer forced to 1 above `w_repro_max`; its value is
  determined entirely by the maturity ogive and the reproductive
  proportion.

Edition 3 keeps all the changes of edition 2 and makes two more:

- The sums over the predation kernel in the encounter rate, the
  predation rate, the predation diffusion and the diet are evaluated by
  fast Fourier transforms on a zero-padded grid (see
  [`convolve_pred_kernel()`](https://sizespectrum.org/mizer/reference/convolve_pred_kernel.md)).
  In earlier editions the convolutions are circular: the part of the
  feeding kernel of a small fish that reaches below the smallest
  resource size wraps around to the largest sizes, so that larvae are
  fed by the largest fish and are counted as eating them.

- A new model uses the second-order numerical scheme in size, as with
  `second_order_w = TRUE`, instead of the first-order scheme of previous
  mizer versions (see
  [`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)).
  The van Leer flux limiter replaces the upwind flux, whose numerical
  diffusion smears out cohorts and can damp real oscillations away.
  Integrals over the size grid are bin averages instead of left-edge
  sums, which moves unnormalised quantities such as biomass and yield by
  around 10% on a typical grid. You can still choose the first-order
  scheme with the `second_order_w` argument of the function that creates
  the model.

## Examples

``` r
# The edition used for new models
defaults_edition()
#> [1] 1
# The edition of an existing model
defaults_edition(NS_params)
#> [1] 1
```
