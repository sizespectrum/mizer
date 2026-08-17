```r
# ── Steady state ──────────────────────────────────────────────────────────────
params <- steady(params)
params <- steadySingleSpecies(params)   # fast starting spectrum
params <- steadyNewton(params)          # direct solve (experimental)

# ── Calibrate to data (re-run steady() after each) ────────────────────────────
params <- calibrateBiomass(params)      # total biomass  → kappa
params <- matchBiomasses(params)        # per-species biomass
params <- calibrateNumber(params)       # same, for `number_observed` instead
params <- matchNumbers(params)          #   of `biomass_observed`
params <- matchGrowth(params)           # growth → h, gamma, ks, k
params <- steady(params)                # re-converge

# ── Reproduction ──────────────────────────────────────────────────────────────
reproduction_level(params) <- 0.25
reproduction_level(params)              # what the model is currently tuned to

# ── Verify ────────────────────────────────────────────────────────────────────
isSteady(params)                        # TRUE if settled within tolerance
summary(params)                         # includes the biomass-drift verdict
plot(getSteadyResidual(params))         # which species and sizes are still moving
plotSpectra(params)
plotGrowthCurves(params)
plotBiomassObservedVsModel(params)
plotYieldObservedVsModel(params)
```
