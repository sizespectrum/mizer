# Balance the active reproduction law while preserving dynamical rates

Adjusts only parameters authorised by a stored, provider-bound tuning
strategy. By default the current effective recruitment is retained, even
when the initial state is not steady. `rdd_target` supplies explicit
targets for a tuning or calibration operation; unselected species retain
their original recruitment and parameters.

## Usage

``` r
balanceRDD(
  params,
  rdd_target = NULL,
  strategy = NULL,
  preserve = NULL,
  ...,
  effort = params@initial_effort,
  t = 0,
  info_level = default_info_level()
)
```

## Arguments

- params:

  A MizerParams object.

- rdd_target:

  Optional scalar, full vector, or named partial vector of recruitment
  targets. `NULL` preserves current effective recruitment. `NA` leaves a
  species unselected.

- strategy:

  Name of a registered reproduction strategy. By default the provider's
  declared default is used.

- preserve:

  Alternative to `strategy`, selecting the corresponding `preserve_...`
  strategy. Built-ins support `"erepro"`, `"reproduction_level"`, and
  their density parameter where applicable.

- ...:

  Law-specific controls described above, or controls accepted by an
  extension's balancing callback.

- effort:

  Fishing effort at the evaluation state.

- t:

  Time at which to evaluate the rates. Stationarity must be declared by
  the provider; it is never inferred by probing other times.

- info_level:

  Controls information and warnings; see
  [`default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.md).

## Value

The verified MizerParams candidate. Errors leave the input unchanged.

## Details

The complete effective rate pipeline is evaluated before and after the
candidate change. Mizer checks unchanged state, declared parameter
bounds and preserved quantities, and every effective rate entering the
equations. Unsupported providers and infeasible controls are errors.
There is no numerical fallback for an unknown reproduction function.

Built-in laws support explicit `erepro`, their density parameter
(`R_max`, `ricker_b` or `sheperd_b`), or `reproduction_level` in `...`.
Supply at most one of these controls. Shepherd also accepts `sheperd_c`.
Ricker and overcompensating Shepherd inversions use `limb = "ascending"`
by default; choose `limb = "descending"` explicitly for the other finite
solution. Built-in defaults preserve their law-specific reproduction
level, except unbounded Shepherd curves, which preserve `sheperd_b`.
`noRDD` adjusts `erepro`, and `constantRDD` adjusts
`constant_reproduction`.

Extension methods may transform targets and delegate with
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html). The outer
transaction verifies the final result once, after the full method chain.
Methods must accept `...` to pass through runtime context. Providing a
method alone does not declare support: the effective provider also needs
an explicit bound tuning strategy.

## See also

[`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md),
[`reproduction_level()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md),
[`setProcess()`](https://sizespectrum.org/mizer/reference/setProcess.md)

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
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md),
[`use_predation_diffusion()`](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)
