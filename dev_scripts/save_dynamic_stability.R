# Computes the long-running calculations for vignettes/dynamic_stability.Rmd
# and saves the results to vignettes/dynamic_stability.rds
#
# The article shows this code with `eval = FALSE` and displays the objects
# saved here, so the two must be kept in step: if you change a calculation in
# the article, change it here and re-run the script.

devtools::load_all(".")
source("dev_scripts/saved_results.R")

params <- NS_params

## The unstable steady state at effort 1.5 ---------------------------------
params_f15 <- findSteadyState(params, effort = 1.5, solver = "newton")

stab <- getStability(params_f15, effort = 1.5)

growth_rate <- function(dt) {
    rho <- getDiscreteStability(params_f15, effort = 1.5, dt = dt)$spectral_radius
    log(rho) / dt
}
growth_rates <- sapply(c(0.1, 0.5, 1), growth_rate)

## The limit cycle the dynamics settle onto --------------------------------
sim_cycle <- projectUntilSettled(params, effort = 1.5, t_max = 200,
                                 t_per = 0.2, method = "tr_bdf2")
conv_cycle <- attr(sim_cycle, "convergence")
# The biomass array carries everything plot() needs, and is far smaller than
# the MizerSim it came from. `plotBiomass(sim)` is `plot(getBiomass(sim))`.
biomass_cycle <- getBiomass(sim_cycle)

attractor_steady <- attr(projectUntilSettled(params, t_max = 100),
                         "convergence")$attractor

## The linearised cycle ----------------------------------------------------
lcs <- getOscillationModeSim(params_f15, amplitude = 0.1)
biomass_lcs <- getBiomass(lcs)

## The bifurcation diagram -------------------------------------------------
scan <- scanModel(params, scan_values = seq(1.0, 1.5, 0.05),
                  set_func = scanEffort(), value_func = getBiomass,
                  species = "Saithe", method = "tr_bdf2",
                  t_max = 2000, amp_rel_tol = 0.01)

## Where the leading eigenvalue crosses zero -------------------------------
efforts <- c(1.00, 1.05, 1.10, 1.15)
max_real <- numeric(length(efforts))
p_e <- params
for (i in seq_along(efforts)) {
    p_e <- findSteadyState(p_e, effort = efforts[i], solver = "newton")
    max_real[i] <- getStability(p_e, effort = efforts[i])$max_real_part
}

## What the default euler step makes of the same model ---------------------
sim_euler <- projectUntilSettled(params, effort = 1.5, t_max = 200,
                                 t_per = 0.2, method = "euler")
conv_euler <- attr(sim_euler, "convergence")

params_nudged <- params_f15
initialN(params_nudged) <- initialN(params_f15) * 1.05
run <- function(method) {
    projectUntilSettled(params_nudged, effort = 1.5, t_max = 200,
                        t_per = 0.2, method = method)
}
sim_nudged_euler <- run("euler")
sim_nudged_trbdf2 <- run("tr_bdf2")
conv_nudged_euler <- attr(sim_nudged_euler, "convergence")
conv_nudged_trbdf2 <- attr(sim_nudged_trbdf2, "convergence")
steady_nudged_euler <- isSteady(finalParams(sim_nudged_euler))

save_with_report(
    list(
        stab = stab,
        growth_rates = growth_rates,
        conv_cycle = conv_cycle,
        biomass_cycle = biomass_cycle,
        attractor_steady = attractor_steady,
        biomass_lcs = biomass_lcs,
        scan = scan,
        efforts = efforts,
        max_real = max_real,
        conv_euler = conv_euler,
        conv_nudged_euler = conv_nudged_euler,
        conv_nudged_trbdf2 = conv_nudged_trbdf2,
        steady_nudged_euler = steady_nudged_euler
    ),
    "vignettes/dynamic_stability.rds"
)
