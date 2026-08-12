```r
# ── Build ─────────────────────────────────────────────────────────────────────
species_params <- read.csv(
    system.file("extdata", "NS_species_params.csv", package = "mizer"))
params <- newMultispeciesParams(species_params, interaction)
params <- newTraitParams()  # or newCommunityParams(), newSingleSpeciesParams()

# ── Inspect ───────────────────────────────────────────────────────────────────
summary(params)
species_params(params)      # given + calculated, one row per species
getInteraction(params)      # the interaction matrix
gear_params(params)         # the fishing gears
resource_params(params)     # the resource scalars

# ── Save / reload ─────────────────────────────────────────────────────────────
params <- setMetadata(params, title = "Celtic Sea model", ...)
saveParams(params, "model.rds")
params <- readParams("model.rds")
saveSim(sim, "sim.rds")
sim <- readSim("sim.rds")
```
