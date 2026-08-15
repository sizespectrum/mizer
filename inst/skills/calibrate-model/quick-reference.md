```r
# ── Steady state ──────────────────────────────────────────────────────────────
params <- steady(params)
params <- steadySingleSpecies(params)   # fast starting spectrum
params <- steadyNewton(params)          # direct solve (experimental)

# ── Calibrate to data (re-run steady() after each) ────────────────────────────
params <- calibrateBiomass(params)      # total biomass  → kappa
params <- matchBiomasses(params)        # per-species biomass
params <- matchGrowth(params)           # growth → h, gamma, ks, k
params <- steady(params)                # re-converge

# ── Reproduction ──────────────────────────────────────────────────────────────
reproduction_level(params) <- 0.25
reproduction_level(params)              # what the model is currently tuned to

# ── Verify ────────────────────────────────────────────────────────────────────
plotSpectra(params)
plotGrowthCurves(params)
plotBiomassObservedVsModel(params)
plotYieldObservedVsModel(params)
```
