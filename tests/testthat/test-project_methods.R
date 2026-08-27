# Initialise --------------------------------------------------------------
# Snapshots recorded with edition 1; lock params creation to edition 1
withr::local_options(mizer_defaults_edition = 1)

# North sea
params <- newMultispeciesParams(NS_species_params_gears_small, inter_small,
                                n = 2/3, p = 0.7, lambda = 2.8 - 2/3,
                                info_level = 0)
no_gear <- dim(params@catchability)[1]
no_sp <- dim(params@catchability)[2]
no_w <- length(params@w)
no_w_full <- length(params@w_full)
sim <- project(params, effort = 1, t_max = 2, dt = 0.5, t_save = 0.5)

# Rescaled
params_r <- params
volume <- 1e-13
params_r@initial_n <- params@initial_n * volume
params_r@initial_n_pp <- params@initial_n_pp * volume
for (res in names(params@initial_n_other)) {
    params_r@initial_n_other[[res]] <- params@initial_n_other[[res]] * volume
}

params_r@species_params$gamma <- params@species_params$gamma / volume
params_r <- setSearchVolume(params_r)
params_r@species_params$R_max <- params_r@species_params$R_max * volume

# Random abundances
set.seed(0)
n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w))) * 1e9
n_full <- abs(rnorm(no_w_full)) * 1e9

# Helper: strip params attribute before snapshot to avoid JSON roundtrip failure
drop_params <- function(x) { attr(x, "params") <- NULL; x }
assign("rate_contribution_at_time", function(params, t, ...) {
    0 * params@mu_b + t
}, envir = globalenv())

params2 <- params
params2@initial_n <- params2@initial_n / 2
params2@initial_n_pp <- params2@initial_n_pp / 2
params2@initial_n_other <- list(test = 1)
params2@initial_effort <- params2@initial_effort / 2


test_that("mizerEReproAndGrowth, mizerERepro and mizerEGrowth follow formulas", {
    encounter <- getEncounter(params, n = n, n_pp = n_full)
    feeding_level <- getFeedingLevel(params, n = n, n_pp = n_full,
                                     encounter = encounter)
    e <- mizerEReproAndGrowth(params, n = n, n_pp = n_full, n_other = list(),
                              t = 0, encounter = encounter,
                              feeding_level = feeding_level)
    expected_e <- sweep((1 - feeding_level) * encounter, 1,
                        params@species_params$alpha, "*",
                        check.margin = FALSE) - params@metab
    expect_equal(e, expected_e, ignore_attr = c("value_name", "units", "class"))

    e_test <- e
    e_test[1, 1] <- -1
    e_repro <- mizerERepro(params, n = n, n_pp = n_full, n_other = list(),
                           t = 0, e = e_test)
    expect_equal(e_repro[1, 1], 0)
    expect_equal(e_repro[-c(1)], (params@psi * pmax(e_test, 0))[-c(1)])

    e_growth <- mizerEGrowth(params, n = n, n_pp = n_full, n_other = list(),
                             t = 0, e_repro = e_repro, e = e_test)
    expect_equal(e_growth, pmax(e_test, 0) - e_repro)
})

test_that("Mortality and encounter contributions receive the current time", {
    e <- getEncounter(params, t = 2)
    m <- getMort(params, t = 2)
    p <- params
    p@other_mort <- list(test = "rate_contribution_at_time")
    expect_equal(getMort(p, t = 2), m + 2, ignore_attr = TRUE)
    p <- params
    p@other_encounter <- list(test = "rate_contribution_at_time")
    expect_equal(getEncounter(p, t = 2), e + 2, ignore_attr = TRUE)
})

test_that("mizerFMortGear, mizerFMort, mizerPredMort and mizerResourceMort follow formulas", {
    effort <- setNames(seq_len(nrow(params@catchability)),
                       dimnames(params@catchability)$gear)
    f_gear <- mizerFMortGear(params, effort)
    expected_f_gear <- params@selectivity
    expected_f_gear[] <- effort * c(params@catchability) * c(params@selectivity)
    expect_equal(f_gear, expected_f_gear)

    f_total <- mizerFMort(params, n = n, n_pp = n_full, n_other = list(),
                          t = 0, effort = effort,
                          e_growth = array(0, dim = dim(params@initial_n)),
                          pred_mort = array(0, dim = dim(params@initial_n)))
    expect_equal(f_total, colSums(f_gear))

    pred_rate <- getPredRate(params, n = n, n_pp = n_full)
    idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
    expect_equal(mizerPredMort(params, n = n, n_pp = n_full, n_other = list(),
                               t = 0, pred_rate = pred_rate),
                 t(params@interaction) %*% pred_rate[, idx_sp, drop = FALSE])
    expect_equal(mizerResourceMort(params, n = n, n_pp = n_full, n_other = list(),
                                   t = 0, pred_rate = pred_rate),
                 as.vector(params@species_params$interaction_resource %*% pred_rate))
})

test_that("mizerRates returns the standard rate list from registered functions", {
    rates_fns <- lapply(params@rates_funcs, get)
    r <- mizerRates(params,
                    n = params@initial_n,
                    n_pp = params@initial_n_pp,
                    n_other = params@initial_n_other,
                    effort = params@initial_effort,
                    t = 0,
                    rates_fns = rates_fns)
    expect_named(r,
                 c("encounter", "feeding_level", "e", "e_repro", "e_growth",
                   "diffusion", "pred_rate", "pred_mort", "f_mort", "mort",
                   "rdi", "rdd", "resource_mort"))
    expect_identical(r$encounter,
                     rates_fns$Encounter(params,
                                         n = params@initial_n,
                                         n_pp = params@initial_n_pp,
                                         n_other = params@initial_n_other,
                                         t = 0))
})

# fft ----
test_that("Test that fft based integrator gives similar result as old code", {
    # Make it harder by working with kernels that need a lot of cutoff
    species_params <- NS_species_params_gears_small
    species_params$pred_kernel_type <- "truncated_lognormal"
    species_params$sigma[2] <- 3
    species_params$beta[1] <- species_params$beta[1] * 100
    species_params$beta[3] <- species_params$beta[3] / 1000
    # and use different egg sizes
    species_params$w_min <- seq(0.001, 1, length.out = no_sp)
    params <- newMultispeciesParams(species_params, inter_small,
                                        no_w = 30, min_w_pp = 1e-12,
                                        info_level = 0)
    # create a second params object that does not use fft
    params2 <- setPredKernel(params, pred_kernel = pred_kernel(params))
    # Test encounter rate integral
    efft <- getEncounter(params, params@initial_n, params@initial_n_pp)
    e <- getEncounter(params2, params@initial_n, params@initial_n_pp)
    # Only check values at fish sizes
    fish <- outer(1:no_sp, 1:no_w, function(i, a) a >= params@w_min_idx[i])
    expect_equal(efft[fish], e[fish], tolerance = 3e-14, ignore_attr = TRUE)
    # Test available energy integral
    prfft <- getPredRate(params, params@initial_n, params@initial_n_pp)
    pr <- getPredRate(params2, params@initial_n, params@initial_n_pp)
    # Due to problem with fft on M1mac, skip this test on CRAN
    skip_on_cran()
    expect_equal(prfft, pr, tolerance = 3e-14, ignore_attr = TRUE)
})
