# Bin average of a power law over geometric bins

Computes the exact average of the power law \\w^d\\ over each bin
\\\[w_j, w\_{j+1}\]\\, i.e. \$\$\overline{w^d}\_j = \frac{1}{\Delta
w_j}\int\_{w_j}^{w\_{j+1}} w^d\\ dw.\$\$

## Usage

``` r
power_law_bin_average(w, dw, d, w_max = Inf)
```

## Arguments

- w:

  Numeric vector of left bin edges \\w_j\\.

- dw:

  Numeric vector of bin widths \\\Delta w_j\\ (same length as `w`).

- d:

  Single numeric exponent of the power law.

- w_max:

  Optional upper cutoff. The power law is taken to be zero above
  `w_max`, with the straddling bin getting the partial average. Defaults
  to `Inf` (no cutoff), in which case the result is identical to the
  uncut formula.

## Value

A numeric vector (same length as `w`) of bin averages of \\w^d\\
(truncated at `w_max` when supplied).

## Details

The integral has a closed form, so the result is exact (not merely
second order): \$\$\overline{w^d}\_j = \frac{w\_{j+1}^{d+1} -
w_j^{d+1}}{(d+1)\\\Delta w_j}, \quad d \neq -1,\$\$
\$\$\overline{w^d}\_j = \frac{\ln(w\_{j+1}/w_j)}{\Delta w_j}, \quad d =
-1.\$\$

This is used by the bin-averaged (second-order) code paths that need the
average of a power-law rate over each bin, for example
[`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md)
and the resource capacity and rate in
[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md).
The grid does not need to be geometric; only the left bin edges `w` and
the bin widths `dw` are used, with \\w\_{j+1} = w_j + \Delta w_j\\.

An optional upper cutoff `w_max` handles a knife-edge truncation of the
power law (for example the resource carrying capacity, which is \\\kappa
w^{-\lambda}\\ below `w_pp_cutoff` and zero above it). The bin
straddling the cutoff then receives the *partial* bin-average — the
power-law average over the part of the bin below `w_max`, divided by the
full bin width — and bins entirely above the cutoff get zero.
