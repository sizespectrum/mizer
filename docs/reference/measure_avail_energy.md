# Measure the available energy in a power-law prey spectrum

Puts the model into the reference state used for deriving the default
`gamma` and `f0`: no fish, and an untruncated resource spectrum given by
the power law \\\kappa w^{-\lambda}\\ across the full `w_full` grid. It
then measures the resulting encounter rate at the largest size and
divides out its power-law dependence on size, so that the returned
number for each species is the coefficient \\A_i\\ in \\E_i(w) = A_i\\
w^{2 + q_i - \lambda}\\.

## Usage

``` r
measure_avail_energy(params, ...)
```

## Arguments

- params:

  A MizerParams object

- ...:

  Unused.

## Value

A vector with one number for each species

## Details

The measurement uses
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md)
rather than
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md).
What is wanted here is a property of the species parameters, not of the
model's dynamics, so it must not go through the extension dispatch chain
or through a rate function that the user has registered with
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md).
An extension that scales the search volume would otherwise fold its own
factor into mizer's `gamma` default, and because that default then
determines the search volume, the factor would be re-applied on every
rebuild (issue \#577). Additive encounter contributions from
`ext_encounter` and `other_encounter` are excluded for the same reason
(issue \#586).

The caller is responsible for having set the `q` and `gamma` species
parameters and the `search_vol` slot that the measurement is to use.

## For extension authors

This function is the hook by which an extension package that changes the
resource structure of the model keeps
[`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md)
and
[`get_f0_default()`](https://sizespectrum.org/mizer/reference/get_f0_default.md)
working. Mizer's own method measures the energy available from the
single resource spectrum in the `initial_n_pp` slot, placed on an
untruncated power law \\\kappa w^{-\lambda}\\ across `w_full` (ignoring
`w_pp_cutoff` so that evaluating at the largest fish size `w[length(w)]`
recovers the scale-free coefficient \\A_i\\ without upper-cutoff
distortion). An extension that replaces that resource by a system of its
own registers a method for its marker class, for example

    measure_avail_energy.mizerMR <- function(params, ...) {
        # Put each of the extension's resources into the untruncated power-law
        # reference state and return the coefficient A_i for each species.
    }

and both defaults then measure against that system without the extension
having to reimplement the rest of the calculation. A method must return
the coefficient \\A_i\\ described above (evaluated against untruncated
power-law resource spectra so that it does not vanish at
`w[length(w)]`), not the encounter rate itself. If your extension
instead calibrates `gamma` against size-truncated resource spectra at a
species-specific reference size (such as `w_inf` or `w_min`), register
methods for
[`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md)
and
[`get_f0_default()`](https://sizespectrum.org/mizer/reference/get_f0_default.md)
directly.

In either case, follow mizer in measuring a property of the species
parameters: leave out any dynamic modulation of the search volume, and
any additive encounter contribution, that your extension applies.
Folding such a factor into the measurement folds it into `gamma`, which
then determines the search volume, so the factor is re-applied on every
rebuild (issues \#577 and \#586).

## See also

Other functions calculating defaults:
[`get_f0_default()`](https://sizespectrum.org/mizer/reference/get_f0_default.md),
[`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md),
[`get_h_default()`](https://sizespectrum.org/mizer/reference/get_h_default.md),
[`get_ks_default()`](https://sizespectrum.org/mizer/reference/get_ks_default.md)
