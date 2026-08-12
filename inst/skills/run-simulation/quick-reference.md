```r
# ── Run ───────────────────────────────────────────────────────────────────────
sim <- project(params, t_max = 20, effort = 1)
sim <- project(params, t_max = 50, dt = 0.01, t_save = 0.5)

# ── Effort ────────────────────────────────────────────────────────────────────
project(params, effort = 1)                        # all gears, constant
project(params, effort = c(Otter = 0.5, Beam = 1)) # per gear, constant
project(params, effort = effort_array)             # time × gear array

gears <- names(initial_effort(params))
years <- 2010:2030
effort_array <- array(1, dim = c(length(years), length(gears)),
                      dimnames = list(time = years, gear = gears))
effort_array[as.character(2020:2030), "Otter"] <- 1.5

# ── Continue / branch a run ───────────────────────────────────────────────────
sim2       <- project(sim, t_max = 10, effort = 2)   # resume from the end
params_end <- finalParams(sim)                       # state at the last step
params_t   <- getParams(sim, time_range = 2010:2015) # averaged over a range

# ── Dynamics studies: avoid numerical diffusion ───────────────────────────────
params <- newMultispeciesParams(sp, second_order_w = TRUE)
second_order_w(params) <- TRUE                       # or on an existing model
sim    <- project(params, t_max = 200, method = "tr_bdf2")

# ── Scenario comparison ───────────────────────────────────────────────────────
sim_low  <- project(params, t_max = 30, effort = 0.5)
sim_high <- project(params, t_max = 30, effort = 1.5)
plotSpectra2(sim_low, sim_high, "F = 0.5", "F = 1.5")
```
