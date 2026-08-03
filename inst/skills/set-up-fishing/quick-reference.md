```r
# ── Gears and selectivity ─────────────────────────────────────────────────────
gp <- gear_params(params)
gp$sel_func <- "sigmoid_length"
gp$l50 <- 25; gp$l25 <- 20
gp$catchability <- 1
gear_params(params) <- gp

# ── Effort ────────────────────────────────────────────────────────────────────
initial_effort(params) <- c(Otter = 0.5, Beam = 1)   # baseline
sim <- project(params, effort = 1)                   # constant during run
sim <- project(params, effort = effort_array)        # time × gear array

# ── Inspect ───────────────────────────────────────────────────────────────────
initial_effort(params)
plotFMort(params)
getFMortGear(params)
```
