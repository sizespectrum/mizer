# Integrating over the size grid

Use this skill when writing anything that integrates over the size grid: a
summary or indicator function, a diagnostic derived from a rate, a new rate
setter with a size-dependent parameter, or an extension that reaches inside the
predation convolution.

mizer has **two quadrature schemes**, selected by the `bin_average` entry of the
`second_order_w` slot. `FALSE` is the default, so code that ignores the flag
looks correct, passes its tests, and is silently wrong by ~10% for every user
who has switched second order on.

## The one rule

Each bin integral is performed in exactly one place. A size-dependent factor is
bin-averaged **where its integral is performed, and nowhere else**.

The theory — why, and where each of mizer's integrals is performed — is in
`vignettes/numerical_details.qmd`, sections *Point values and bin averages* and
*The `second_order_w` switch*. Read those rather than re-deriving; the tables
there are the authoritative inventory and are meant to be kept current.

## Decide which case you are in

### Case 1 — a plain integral against the abundance, ∫ K(w) N(w) dw

Discretise as `sum_j Kbar_j * N_j * dw_j`. Only `K` is approximated:

```r
K <- bin_average_summary_weight(K, params)   # gated on the flag
drop(n %*% (K * params@dw))
```

- **Gate it, don't hard-code it.** `bin_average_summary_weight()` returns `K`
  untouched on the default path, so the old numbers stay byte-identical.
- **Never bin-average `N` or `dw`.** `N_j` is already a cell average and `dw_j`
  is exact.
- **Average the product, not the factors.** SSB averages `psi * w`; yield
  averages `F * w`. Averaging separately is a different (wrong) number.
- **If `K` is an exact power law `w^a`,** use `power_law_bin_average(w, dw, a)`
  instead of the trapezoid — it is exact, not merely second order.
- If the result is size-resolved, tag it: `ArraySpeciesBySize(..., representation
  = "average")` for a bin average, `"point"` for a boundary quantity. The tag
  drives the half-bin plotting shift.

### Case 2 — a quantity built from rates mizer already computes

Call `getEncounter()`, `getFeedingLevel()`, `getPredRate()`, `getEGrowth()`, …
and do not rebuild them. The rate functions already carry the right quadrature
for whichever scheme the model is in. Re-deriving a rate is how the two known
bugs in this area were introduced.

### Case 3 — you need to go inside the encounter or predation convolution

Use **`encounter_kernel(params)`**, not `getPredKernel(params)`, and pair it
with the **plain point weight** `params@w * params@dw`.

Under `bin_average`, `setPredKernel()` builds `ft_pred_kernel_e` from the kernel
integrated over the prey bin and divides by `beta - 1` precisely so that the
`w * dw` supplied by the prey vector cancels. That `w * dw` is a normalisation,
not a first-order quadrature weight — bin-averaging it applies the prey-bin
integral twice.

### Case 4 — a new rate setter with a size-dependent parameter

Follow `setExtMort()` / `setExtDiffusion()` / `setResource()`: gate on the flag
and use `power_law_bin_average()` for power laws, or a composite midpoint rule
(as `setFishing()` does for selectivity) for anything else. Do the integral once,
at setup, so the projection cost is unchanged.

## Traps

### Double-counting is a uniform factor, so normalised outputs hide it

On a geometric grid `bin_average_weight(w) / w` is exactly `(1 + beta) / 2`
(1.0967 for `NS_params`). Applying the prey-bin quadrature twice therefore
scales the result by a constant, which **cancels in any proportion or ratio**.
`getDiet(proportion = FALSE)` was 9.7% too large for a long time while the
default `proportion = TRUE` stayed correct (#474). If a consistency ratio comes
out as a constant, read off its value: `(1 + beta) / 2` means the quadrature was
applied twice, `2 / (1 + beta)` means it is missing.

### `getPredKernel()` is not the kernel the encounter uses

It returns the kernel point-sampled on the grid — right for plotting, and the
form you supply a custom kernel in, but not the bin-integrated coefficients the
convolution consumes. Pairing it with `getEncounter()` in a numerator/denominator
is what made `getTrophicLevel()` wrong (#474).

### Growth-type rates are never bin-averaged

`g`, `e`, the encounter rate and the feeding level are point values at `w_j`
under **both** settings — they are boundary velocities. What improves them when
the flag is on is the encounter integral behind them, not any averaging of the
rate itself. Bin-averaging them is an error, not an upgrade.

### Testing only the default path proves nothing

Both #474 bugs were invisible with `bin_average = FALSE`. Every new integral
needs a test with the flag on.

## Verifying a change

1. **Default path unchanged.** Assert byte-identity (or an existing snapshot)
   with `bin_average = FALSE`. Any movement there is a regression.
2. **Flag on: assert the identity your quantity should satisfy.** Anything that
   decomposes a rate must reassemble into it:

   ```r
   params <- NS_params_small
   second_order_w(params) <- c(bin_average = TRUE)
   total <- rowSums(getDiet(params, proportion = FALSE), dims = 2)
   ratio <- total / (getEncounter(params) * (1 - getFeedingLevel(params)))
   range(ratio[initialN(params) > 0])   # 1 1
   ```

   FFT convolution is circular, so allow ~1e-4 when comparing a direct sum
   against a rate function; within one code path the agreement is exact.
3. **Convergence.** The gap between the two schemes should shrink under grid
   refinement — see the "second-order biomass converges to default" test for the
   pattern.

New tests go in `tests/testthat/test-second_order_summary.R`, using the
`NS_params_small` fixture and toggling with
`second_order_w(p) <- c(bin_average = TRUE)`.

## Helpers

| Helper | Use |
|---|---|
| `bin_average_summary_weight(K, params)` | trapezoidal bin average, gated on the flag — the default entry point |
| `bin_average_weight(K)` | ungated trapezoid; averages along the last dimension of an array |
| `power_law_bin_average(w, dw, a, w_max)` | exact bin average of `w^a`, with optional cutoff |
| `encounter_kernel(params)` | the kernel `mizerEncounter()` actually uses, under either scheme |
| `bin_midpoints(params)` | geometric bin centres, for plotting bin averages |

Setting `second_order_w(params) <- c(bin_average = TRUE)` re-runs `setParams()`,
because every array in the inventory table is precomputed.
