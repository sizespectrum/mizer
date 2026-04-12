# Upgrade MizerParams / MizerSim data objects

Use this skill after any change to the `MizerParams` or `MizerSim` S4 class
definition — new slots, removed slots, changed slot types, or new entries in
named lists such as `@rates_funcs`.

## Steps

### 1. Update `upgradeParams()` in `R/upgrade.R`

Add code that patches old objects to the new state. Place it just before the
final block:

```r
params@mizer_version <- packageVersion("mizer")
params <- validParams(params, info_level = 0)
```

Example — adding a new `rates_funcs` entry:

```r
# Add Diffusion rate function if missing (added in 2.5.4.9122)
if (is.null(params@rates_funcs[["Diffusion"]])) {
    params@rates_funcs[["Diffusion"]] <- "mizerDiffusion"
}
```

### 2. Bump the version threshold in `needs_upgrading()` in `R/upgrade.R`

Change the comparison to the **new** version so that existing objects stored
with the old version are detected as needing an upgrade:

```r
!.hasSlot(params, "mizer_version") ||
    params@mizer_version < "2.5.4.9122"   # ← new version
```

### 3. Bump the version in `DESCRIPTION`

```
Version: 2.5.4.9122
```

### 4. Reload, upgrade, and save the data objects

```r
devtools::load_all()
NS_params <- upgradeParams(NS_params)
NS_sim    <- upgradeSim(NS_sim)
save(NS_params, file = "data/NS_params.rda", compress = "xz")
save(NS_sim,    file = "data/NS_sim.rda",    compress = "xz")
```

### 5. Reload again before running tests

```r
devtools::load_all()
devtools::test()
```

## Why this matters

Tests load `NS_params` and `NS_sim` from the `data/` directory. If those
objects predate the class change they will be missing new slots or list
entries, causing errors that look like rate-function or projection failures
rather than the real cause (stale stored objects).
