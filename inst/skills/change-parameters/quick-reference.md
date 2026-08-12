```r
# ── Species parameters (preferred: triggers recalculation) ────────────────────
species_params(params)$beta <- 150
species_params(params)                     # view all (given + calculated)
calculated_species_params(params)          # only what mizer derived

# ── Size-dependent rate arrays ────────────────────────────────────────────────
params <- setMetabolicRate(params)         # recompute from species params
metab(params) <- my_array                  # direct setter; freezes the array
params <- setSearchVolume(params, search_vol = my_array)  # same, via set...()
params <- setSearchVolume(params, reset = TRUE)  # unfreeze: recompute again
params <- setParams(params)                # rebuild ALL rate arrays

# ── Fishing ───────────────────────────────────────────────────────────────────
gear_params(params) <- gp
params <- setFishing(params, initial_effort = c(Otter = 1))

# ── Resource (scalars rebuild the arrays, like species parameters) ────────────
resource_params(params)$kappa <- 1e11      # rescales the carrying capacity
resource_capacity(params) <- my_capacity   # bespoke array; freezes it
params <- setResource(params, reset = TRUE)  # unfreeze: recompute from scalars

# ── Interaction ───────────────────────────────────────────────────────────────
params <- setInteraction(params, inter)

# ── Rate functions (level 3: replace how a rate is computed) ──────────────────
seasonalMort <- function(params, t, ...) {
    mizerMort(params, t = t, ...) * (1 + 0.3 * sin(2 * pi * t))   # t in years
}
params <- setRateFunction(params, "Mort", "seasonalMort")  # now time-dependent
other_params(params)$my_param <- 42                        # params for your function
```
