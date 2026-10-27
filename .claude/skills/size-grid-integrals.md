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

**Call `sizeIntegral()`.** It is the one place this integral is written:

```r
sizeIntegral(params, weighting = K, min_w = 10, max_w = 5000)
```

It applies the size-range mask, gates the bin-averaging on the flag, multiplies
by `dw`, contracts over the size axis and wraps the result. `getBiomass()`,
`getN()`, `getSSB()`, `getYield()`, `getYieldGear()` and
`getProportionOfLargeFish()` are all implemented with it; add the next one the
same way rather than writing the sum again. Do not pass `params@dw` in the
weighting factor and do not call `bin_average_weight()` before handing the
weighting factor over — both are done inside.

- **Average the product, not the factors.** Pass the whole product as `weighting`:
  SSB uses `psi * w`; yield uses `F * w`. Averaging separately is a different
  (wrong) number. The size-range mask counts as one of the factors — it is
  multiplied in *before* the averaging, which is what makes the bin straddling
  the boundary contribute partially.
- **Never bin-average `N` or `dw`.** `N_j` is already a cell average and `dw_j`
  is exact.
- **If `K` is an exact power law `w^a`,** and you are writing a rate setter
  rather than a summary, use `power_law_bin_average(w, dw, a)` instead of the
  trapezoid — it is exact, not merely second order.
- If the result is size-resolved it is not this case, so `sizeIntegral()` does
  not apply; tag it yourself: `ArraySpeciesBySize(..., representation =
  "average")` for a bin average, `"point"` for a boundary quantity. The tag
  drives the half-bin plotting shift.

The contraction inside `sizeIntegral()` is a matrix multiplication, chosen so
that it reproduces mizer's historical `n %*% (K * dw)` to the last bit. Changing
it to `rowSums()` (long-double accumulation) moves results by ~1e-16 and breaks
that guarantee.

### Case 2 — a quantity built from rates mizer already computes

Call `getEncounter()`, `getFeedingLevel()`, `getPredRate()`, `getEGrowth()`, …
and do not rebuild them. The rate functions already carry the right quadrature
for whichever scheme the model is in. Re-deriving a rate is how the two known
bugs in this area were introduced.

### Case 3 — you need to go inside the encounter or predation convolution

Use **`encounter_kernel(params)`**, not `pred_kernel(params)`, and pair it
with the **plain point weight** `params@w * params@dw`.

Under `bin_average`, `setPredKernel()` builds `ft_pred_kernel_e` from the kernel
integrated over the prey bin and divides by `beta - 1` precisely so that the
`w * dw` supplied by the prey vector cancels. That `w * dw` is a normalisation,
not a first-order quadrature weight — bin-averaging it applies the prey-bin
integral twice.

The predation convolution is the transpose, and there the opposite holds for the
predator-side factor. `ft_pred_kernel_p` cancels only the predator `dw`; the
predator's `(1 - f) * search_vol` is a smooth weight of that integral, so the
predator vector is
`bin_average_weight((1 - feeding_level) * search_vol, params) * n * dw`, as in
`mizerPredRate()` (#619). The kernel and that factor are averaged separately,
unlike the single weight of Case 1, because the feeding level changes every
step and cannot go into the precomputed kernel. The covariance this ignores is
O(Δx²), so it keeps second order.

### Case 4 — a new rate setter with a size-dependent parameter

Follow `setExtMort()` / `setExtDiffusion()` / `setResource()`: gate on the flag
and use `power_law_bin_average()` for power laws, or a composite midpoint rule
(as `setFishing()` does for selectivity) for anything else. Do the integral once,
at setup, so the projection cost is unchanged.

## Traps

### Double-counting is a uniform factor, so normalised outputs hide it

On a geometric grid `trapezoidal_bin_average(w) / w` is exactly `(1 + beta) / 2`
(1.0967 for `NS_params`). Applying the prey-bin quadrature twice therefore
scales the result by a constant, which **cancels in any proportion or ratio**.
`getDiet(proportion = FALSE)` was 9.7% too large for a long time while the
default `proportion = TRUE` stayed correct (#474). If a consistency ratio comes
out as a constant, read off its value: `(1 + beta) / 2` means the quadrature was
applied twice, `2 / (1 + beta)` means it is missing.

### `pred_kernel()` is not the kernel the encounter uses

It returns the kernel point-sampled on the grid — right for plotting, and the
form you supply a custom kernel in, but not the bin-integrated coefficients the
convolution consumes. Pairing it with `getEncounter()` in a numerator/denominator
is what made `getTrophicLevel()` wrong (#474).

### Growth-type rates are stored as point values — but may be integrands elsewhere

`g`, `e`, the encounter rate, the feeding level and the search volume are
stored as point values at `w_j` under **both** settings. `g` is a boundary
velocity, and `E(w) = gamma(w) * ∫ ...` has `gamma` outside its integral. What
improves them when the flag is on is the encounter integral behind them, not any
averaging of the rate itself. Bin-averaging the stored rates is an error, not an
upgrade.

That does not make them point values in **every** integral they enter. A
factor's role is set by the integral it appears in, relative to that
integral's variable of integration:

| Point value in | Integrand in | Averaged where |
|---|---|---|
| growth (`e`) | reproduction `∫ psi e N dw` | `mizerRDI()`, as the product `psi * e` |
| encounter, growth (`f`, `search_vol`) | predation `∫ phi (1 - f) gamma N dw`, over **predator** size | `mizerPredRate()`, as the product `(1 - f) * search_vol` |
| encounter, growth (`f`, `search_vol`) | predation-diffusion coefficient `(1 - f) gamma alpha^2 I_d`, a predator-bin average | `mizerDiffusion()`, as the product `(1 - f) * search_vol` |

Average the product with `bin_average_weight()` at that use-site, and leave the
stored rate alone. #619 was this trap: predation was derived "in the same way"
as the encounter, `(1 - f) * gamma * N * dw` was treated as "the predator vector"
by analogy with the prey vector, and only the kernel was integrated over the
predator bin. Predation and resource mortality were first order for it. The
diffusion coefficient had the same flaw from the other side: `setPredKernel()`
folds the kernel so that `I_d` comes out averaged over the predator bin, but the
factor in front of it was still point-sampled. Before deciding a factor is a
point value, write the integral out and check whether the factor depends on the
integration variable or sits in front of an integral that is being averaged.

### Testing only the default path proves nothing

Both #474 bugs were invisible with `bin_average = FALSE`. Every new integral
needs a test with the flag on.

### Testing the parts against the derivation proves nothing about the order

The bin-integrated kernel weights are tested against their analytic values,
and those tests passed throughout #619, because the kernel was right and the
error was in a factor outside it. A test written from the derivation can only
confirm that the code matches the derivation. Only a grid-refinement test of the
rate itself checks that the derivation is complete.

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
4. **Order.** For a rate, measure the order directly. Prescribe smooth
   cell-averaged densities, evaluate the rate at a fixed size on grids an octave
   apart (bin averages at the bin centre `sqrt(w * (w + dw))`, point values at
   the node), and take `log2` of the ratio of successive differences: about 1 is
   first order, about 2 second. `convergence_params()`, `at_bin_centre()` and
   `convergence_order()` in `helper.R` do this. The #619 tests in
   `test-project_methods.R` and `test-diffusion.R` are the pattern, and cost
   about a second each on 80/160/320 bins.

New tests go in `tests/testthat/test-second_order_summary.R`, using the
`NS_params_small` fixture and toggling with
`second_order_w(p) <- c(bin_average = TRUE)`. Tests of `sizeIntegral()` and its
helpers themselves belong in `test-sizeIntegral.R`, per the ordinary rule that a
test lives in the file named after the R file that defines the function.

## Helpers

| Helper | Exported? | Use |
|---|---|---|
| `sizeIntegral(object, weighting, ...)` | yes | the whole integral ∫ N K dw — the default entry point |
| `bin_average_weight(K, params)` | yes | trapezoidal bin average, gated on the flag, for a weight you are not integrating |
| `trapezoidal_bin_average(K)` | no | ungated trapezoid; averages along the last dimension of an array |
| `power_law_bin_average(w, dw, a, w_max)` | no | exact bin average of `w^a`, with optional cutoff |
| `encounter_kernel(params)` | yes | the kernel `mizerEncounter()` actually uses, under either scheme |
| `bin_midpoints(params)` | no | geometric bin centres, for plotting bin averages |

`sizeIntegral()`, `bin_average_weight()` and `encounter_kernel()` are exported
(badged experimental, to match `second_order_w()`) so that extension authors and
users writing their own indicators can reach them; the user-facing guidance lives
in `inst/skills/analyse-and-plot/SKILL.md` and `inst/skills/extend-mizer/SKILL.md`. Keep
those two in step with any change here. The other three are internal — prefer
`ArraySpeciesBySize(..., representation = "average")` over telling a user to call
`bin_midpoints()` by hand.

`sizeIntegral()`'s own helpers (`size_dim_labels()`, `merge_dim_labels()`,
`broadcast_dims()`, `dim_extents()`, `collect_dimnames()`, all in
`R/sizeIntegral.R`) identify each dimension of the abundance and the weight by a
label — `"time"`, `"gear"`, `"sp"`, `"w"` — taken from `names(dimnames())`, and
merge the two label sets. That is what lets a `time x gear x sp x w` weight line
its time dimension up with the abundance's instead of multiplying it out. A new
array shape that needs supporting is a change to `size_dim_labels()`, not a new
branch in `sizeIntegral()`.

Setting `second_order_w(params) <- c(bin_average = TRUE)` re-runs `setParams()`,
because every array in the inventory table is precomputed.
