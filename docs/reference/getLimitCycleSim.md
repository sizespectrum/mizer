# Construct a MizerSim of the linearised limit cycle

**\[experimental\]** Using the leading complex eigenvector from
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md),
constructs a
[MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
object covering one period of the limit cycle in the linear
approximation. The result can be inspected with all standard mizer
plotting functions (e.g.
[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md),
[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md)).

## Usage

``` r
getLimitCycleSim(x, amplitude = 0.1, t_save = 0.1, ...)
```

## Arguments

- x:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object at a steady state, typically the output of
  [`steadyNewton()`](https://sizespectrum.org/mizer/reference/steadyNewton.md),
  or the list returned by
  [`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md).

- amplitude:

  Maximum relative perturbation \\\max_w \|\delta N(t,w)\|/N^\*(w)\\
  across the limit cycle. Default `0.1`.

- t_save:

  The time interval between saved time steps in the returned
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md).
  Defaults to `0.1`.

- ...:

  Additional arguments forwarded to
  [`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
  when `x` is a `MizerParams` object.

## Value

A [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
object whose time axis spans one period \\\[0, T\]\\ of the linearised
limit cycle.

## Details

### Mathematical background

Near a Hopf bifurcation the dominant eigenvalue of the linearised
one-step map is a complex conjugate pair \\\lambda = r e^{i\theta}\\
with \\r \approx 1\\ and angular frequency \\\theta = 2\pi/T\\ per time
step. The linearised perturbation is \$\$\delta N(t) =
A\\\operatorname{Re}\[e^{i\theta t}\\\mathbf{v}\],\$\$ where
\\\mathbf{v}\\ is the leading complex eigenvector (normalised so
\\\max_w \|\mathbf{v}(w)\| = 1\\) and \\A\\ is chosen so that the
maximum *relative* perturbation \\\max\_{t,w}\|\delta N(t,w)\|/N^\*(w) =
\\ `amplitude`. The full state at each time step is \$\$N(t) =
\max(N^\* + \delta N(t),\\ 0)\$\$ (clipping prevents negative abundances
at large amplitudes).

The resource is placed at its quasi-static semichemostat equilibrium
\\n\_{pp}^\*(N(t))\\ at each step, consistent with the reduced Jacobian
used by
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md).

The returned
[MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
has times running from 0 to \\T\\ (the period in time steps, typically
years).

## See also

[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md),
[`steadyNewton()`](https://sizespectrum.org/mizer/reference/steadyNewton.md)
