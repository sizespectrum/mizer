# Construct a MizerSim of the leading oscillatory mode

**\[experimental\]** Using the leading complex eigenvector from
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md),
constructs a
[MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
object covering one period of that oscillation in the linear
approximation. The result can be inspected with all standard mizer
plotting functions (e.g.
[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md),
[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md)).

## Usage

``` r
getOscillationModeSim(x, amplitude = 0.1, t_save = 0.1, ...)
```

## Arguments

- x:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object at a steady state, typically the output of
  [`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md),
  or the list returned by
  [`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md).

- amplitude:

  Largest relative swing in species biomass across the cycle, \\\max_i
  \max_t \|B_i(t) - B_i^\ast\| / B_i^\ast\\. Default `0.1`, meaning the
  most strongly oscillating species departs 10 % from its steady
  biomass.

- t_save:

  The time interval between saved time steps in the returned
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md).
  Defaults to `0.1`. The final interval is shorter when `t_save` does
  not divide the period, so that the cycle closes.

- ...:

  Additional arguments forwarded to
  [`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
  when `x` is a `MizerParams` object.

## Value

A [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
object whose time axis spans one period \\\[0, T\]\\ of the linearised
oscillatory mode.

## Details

The object shows the *shape* of the mode — which species swing, how far,
and in what phase relative to each other and to the resource. Whether
the model actually settles onto this oscillation is a separate question,
answered by the real part of the eigenvalue: it is a limit cycle only
where that real part is zero, at a Hopf bifurcation.

### Mathematical background

An oscillatory mode is a complex-conjugate pair of eigenvalues \\\lambda
= \sigma \pm i\omega\\ of the Jacobian, with period \\T =
2\pi/\|\omega\|\\.
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
returns the pair with the largest \\\sigma\\ as
`leading_oscillatory_eigenvalue` and its eigenvector as
`leading_oscillatory_eigenvector`. The linearised perturbation of the
full state \\x = (N, n\_{pp})\\ is \$\$\delta x(t) =
A\\\operatorname{Re}\[e^{i\omega t}\\\mathbf{v}\],\$\$ where
\\\mathbf{v}\\ is that eigenvector and \\A\\ is chosen so that the
largest *relative swing in species biomass* equals `amplitude`. Biomass
is a linear functional of the abundance, so \$\$B_i(t) = B_i^\* +
A\\\operatorname{Re}\[e^{i\omega t} c_i\], \qquad c_i = \int v_i(w)\\ w
\\ dw,\$\$ and species \\i\\ departs from its steady biomass by at most
\\A\|c_i\|\\. \\A\\ is set so that \\\max_i A\|c_i\|/B_i^\*\\ is
`amplitude`: no species' biomass moves further than that fraction from
its steady value, and the one that oscillates hardest moves exactly that
far. The integral uses
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md),
so it follows the model's own quadrature scheme and agrees with
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
bin for bin.

The cap is on the species that swings hardest rather than on the
community total, because species oscillating out of phase cancel in the
total: a modest total swing can be produced by wild swings in the
individual species.

The state at each time is \$\$x(t) = \max(x^\* + \delta x(t),\\ 0).\$\$
Because the cap is on biomass, an individual size class can still be
driven negative while the biomass it belongs to moves only a little — a
cohort trough is a much larger relative excursion than the biomass
integral over it. That clipping is reported when it happens, and means
the picture is no longer the linear mode.

The fish and resource blocks of \\\mathbf{v}\\ carry a single common
normalisation, so the same \\A\\ drives both and the resource oscillates
with the amplitude and phase *the mode gives it* — generally neither in
step with the fish nor slaved to them. `amplitude` is set on the fish
biomass, so how far the resource moves is a property of the mode rather
than something you choose.

The growth of the mode is deliberately dropped: \\e^{\sigma t}\\ is
omitted so that the oscillation closes after one period. That is exact
only at a Hopf bifurcation, where \\\sigma = 0\\; away from it the
returned object shows the *shape* of the oscillation, not its envelope.
\\\sigma\\ is recorded in the result's `sim_params` as `growth_rate`,
and it is the number to look at before calling what you are seeing a
cycle.

The returned
[MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
has times running from 0 to exactly \\T\\ (the period, in years). The
saved times are spaced `t_save` apart, except for the last interval,
which is shortened when `t_save` does not divide \\T\\. Ending exactly
at \\T\\ is what makes the cycle close: the phase factor \\e^{i\omega
T}\\ is 1, so the final state is the first state again.

## See also

[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md),
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
