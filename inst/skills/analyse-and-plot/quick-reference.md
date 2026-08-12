```r
# ── Accessing raw arrays ───────────────────────────────────────────────────────
N(sim)                  # time × species × size
NResource(sim)          # time × size
finalN(sim)             # species × size  (last time step)
finalNResource(sim)     # size            (last time step)

# ── Species biomass / abundance / yield (time × species) ──────────────────────
getBiomass(sim)         # total biomass
getSSB(sim)             # spawning stock biomass
getN(sim)               # total abundance (numbers)
getYield(sim)           # catch in weight
getYieldGear(sim)       # catch by gear (time × gear × species)

# ── Rates at size (time × species × size) ─────────────────────────────────────
getFeedingLevel(sim)    # satiation (0 = starving, 1 = full)
getPredMort(sim)        # predation mortality
getFMort(sim)           # fishing mortality
getFMortGear(sim)       # fishing mortality by gear (time × gear × species × size)

# ── Diet and trophic (species × size × …) ────────────────────────────────────
getDiet(params)                  # proportion of diet from each prey
getTrophicLevel(params)          # trophic level at size (species × size)
getTrophicLevelBySpecies(params) # mean trophic level (species)

# ── Community indicators (time series) ────────────────────────────────────────
getProportionOfLargeFish(sim, threshold_w = 100)
getMeanWeight(sim)
getMeanMaxWeight(sim)
getCommunitySlope(sim)          # returns data.frame with slope, intercept, R²

# ── Dedicated plot functions ──────────────────────────────────────────────────
# Each plot*() is a shortcut for plot() on the matching get*() array, and each has
# an interactive plotly*() twin (plotlyBiomass(), plotlySpectra(), …).
plot(sim)               # 5-panel summary
plotBiomass(sim)        # biomass vs time
plotYield(sim)          # yield vs time
plotYieldGear(sim)      # yield vs time, faceted by gear
plotSpectra(sim)        # abundance spectra vs size (+ resource & background)
plotFeedingLevel(sim)   # feeding level vs size
plotPredMort(sim)       # predation mortality vs size
plotFMort(sim)          # fishing mortality vs size
plotGrowthCurves(sim)   # size vs age
plotDiet(params, species = "Cod")  # diet composition vs size
plotCDF(sim)            # cumulative biomass/abundance over size

# ── Plot any array directly, plus combine / compare tools ─────────────────────
plot(getResourceMort(params))   # any get*() array plots directly
p <- plot(getBiomass(sim), species = "Cod")
addPlot(p, getBiomass(sim), species = "Herring", linetype = "dashed")  # add lines
plot2(getFMort(params), getFMort(params2), "Before", "After")  # compare arrays
plotRelative(getEGrowth(params), getEGrowth(params2))          # relative diff
plotHover(getBiomass(sim))      # interactive (hover) version of an array plot
animate(sim)                    # animate spectra through time

# ── Compare two simulations or models ─────────────────────────────────────────
plotSpectra2(params, params2, "Before", "After")
plotSpectraRelative(params, params2)      # relative difference of spectra
plotCDF2(sim, sim2, "Unfished", "Fished")
```
