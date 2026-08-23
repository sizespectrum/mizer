# Computes the long-running calculations for vignettes/discontinuous_rates.Rmd
# and saves the results to vignettes/discontinuous_rates.rds

devtools::load_all(".")
source("dev_scripts/saved_results.R")

params <- NS_params
initial_effort(params)["Otter"] <- 1.5
ssb0 <- getSSB(params)[["Cod"]]

other_params(params)$stock     <- "Cod"
other_params(params)$gear      <- "Otter"
other_params(params)$b_lim     <- 0.5 * ssb0
other_params(params)$b_trigger <- 0.7 * ssb0

ssb_of <- function(params, n, species) {
    initialN(params) <- n
    getSSB(params)[[species]]
}

hardHCR <- function(params, n, n_pp, n_other, t, effort, ...) {
    op <- other_params(params)
    if (ssb_of(params, n, op$stock) < op$b_lim) {
        effort[op$gear] <- 0
    }
    mizerFMort(params, n, n_pp, n_other, t = t, effort = effort, ...)
}

rampHCR <- function(params, n, n_pp, n_other, t, effort, ...) {
    op <- other_params(params)
    ssb <- ssb_of(params, n, op$stock)
    frac <- (ssb - op$b_lim) / (op$b_trigger - op$b_lim)
    effort[op$gear] <- effort[op$gear] * min(1, max(0, frac))
    mizerFMort(params, n, n_pp, n_other, t = t, effort = effort, ...)
}

params_hard <- setRateFunction(params, "FMort", "hardHCR")
params_ramp <- setRateFunction(params, "FMort", "rampHCR")

ssb_series <- function(p, dt, t_max = 30, method = "tr_bdf2") {
    sim <- project(p, dt = dt, t_max = t_max, t_save = dt, method = method,
                   progress_bar = FALSE)
    ssb <- getSSB(sim)[, "Cod"]
    data.frame(t = as.numeric(names(ssb)), ssb = ssb / ssb0)
}

hard_coarse <- ssb_series(params_hard, dt = 0.1)
hard_fine   <- ssb_series(params_hard, dt = 0.0125)
ramp_coarse <- ssb_series(params_ramp, dt = 0.1)
ramp_fine   <- ssb_series(params_ramp, dt = 0.0125)

tally <- new.env()

tallyHCR <- function(params, n, n_pp, n_other, t, effort, ...) {
    op <- other_params(params)
    closed <- ssb_of(params, n, op$stock) < op$b_lim
    tally$decisions <- c(tally$decisions, closed)
    if (closed) effort[op$gear] <- 0
    mizerFMort(params, n, n_pp, n_other, t = t, effort = effort, ...)
}
params_tally <- setRateFunction(params, "FMort", "tallyHCR")

switching <- function(dt, t_max = 30) {
    tally$decisions <- logical(0)
    project(params_tally, dt = dt, t_max = t_max, t_save = t_max,
            method = "tr_bdf2", progress_bar = FALSE)
    d <- tail(tally$decisions, length(tally$decisions) %/% 2)
    c(closed_fraction   = mean(d),
      switches_per_year = sum(diff(d) != 0) / (t_max / 2))
}

dts_sw <- c(0.1, 0.05, 0.0125)
sw_table <- cbind(dt = dts_sw, as.data.frame(t(sapply(dts_sw, switching))))

amplitude <- function(p, dt, method) {
    s <- ssb_series(p, dt, method = method)
    tail_s <- s$ssb[s$t > max(s$t) * 2 / 3]
    max(tail_s) - min(tail_s)
}

dts <- c(0.1, 0.05, 0.025, 0.0125)
amp_table <- data.frame(
    dt          = dts,
    euler       = sapply(dts, amplitude, p = params_hard, method = "euler"),
    tr_bdf2     = sapply(dts, amplitude, p = params_hard, method = "tr_bdf2"),
    ramp_tr_bdf2 = sapply(dts, amplitude, p = params_ramp, method = "tr_bdf2")
)

settle <- function(p) {
    sim <- project(p, dt = 0.01, t_max = 40, t_save = 40, progress_bar = FALSE)
    initialN(p) <- finalN(sim)
    initialNResource(p) <- finalNResource(sim)
    p
}

hard_settled <- settle(params_hard)
ramp_settled <- settle(params_ramp)

ramp_ss <- findSteadyState(ramp_settled, solver = "newton")
hard_ss <- suppressWarnings(findSteadyState(hard_settled, solver = "newton"))

ratio_ssb <- c(ramp = getSSB(ramp_ss)[["Cod"]] / other_params(params)$b_lim,
               hard = getSSB(hard_ss)[["Cod"]] / other_params(params)$b_lim)

stab_hard <- suppressWarnings(sapply(c(1e-3, 1e-4, 1e-5), function(h) {
    getStability(hard_ss, h = h)$max_real_part
}))

stab_ramp <- suppressWarnings(sapply(c(1e-3, 1e-4, 1e-5), function(h) {
    getStability(ramp_ss, h = h)$max_real_part
}))

st <- suppressWarnings(getStability(hard_settled))
stab_hard_settled <- c(max_real_part = st$max_real_part, stable = st$stable)

res <- list(
    hard_coarse = hard_coarse,
    hard_fine = hard_fine,
    ramp_coarse = ramp_coarse,
    ramp_fine = ramp_fine,
    sw_table = sw_table,
    amp_table = amp_table,
    hard_settled = hard_settled,
    ramp_settled = ramp_settled,
    ramp_ss = ramp_ss,
    hard_ss = hard_ss,
    ratio_ssb = ratio_ssb,
    stab_hard = stab_hard,
    stab_ramp = stab_ramp,
    st = st,
    stab_hard_settled = stab_hard_settled
)

save_with_report(res, "vignettes/discontinuous_rates.rds")
