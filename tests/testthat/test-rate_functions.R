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

params2 <- params
params2@initial_n <- params2@initial_n / 2
params2@initial_n_pp <- params2@initial_n_pp / 2
params2@initial_n_other <- list(test = 1)
params2@initial_effort <- params2@initial_effort / 2


test_that("getRates works", {
    r <- getRates(params)
    expect_identical(names(r),
                     c("encounter", "feeding_level", "e", "e_repro",
                       "e_growth", "diffusion", "pred_rate", "pred_mort",
                       "f_mort", "mort", "rdi", "rdd", "resource_mort"))
    # test that the optional parameters take the correct defaults
    expect_identical(r,
                     getRates(params, n = params@initial_n,
                              n_pp = params@initial_n_pp,
                              n_other  = params@initial_n_other,
                              effort = params@initial_effort,
                              t = 0))
    # test that getRates actually uses its optional arguments
    expect_identical(getRates(params2),
                     getRates(params, n = params2@initial_n,
                              n_pp = params2@initial_n_pp,
                              n_other  = params2@initial_n_other,
                              effort = params2@initial_effort,
                              t = 0))
})

# getEncounter --------------------------------------------------------------

test_that("getEncounter returns with correct dimnames", {
    enc <- getEncounter(params)
    expect_identical(dimnames(enc),
                     dimnames(params@initial_n))
})

test_that("getEncounter is independent of volume", {
    enc <- getEncounter(params)
    enc_r <- getEncounter(params_r)
    expect_equal(enc, enc_r, ignore_attr = "params")
})

test_that("External encounter is included", {
    enc <- getEncounter(params)
    # add something of the right dimension
    extra_enc <- params@mu_b
    ext_encounter(params) <- ext_encounter(params) + extra_enc
    expect_equal(getEncounter(params), enc + extra_enc, ignore_attr = TRUE)
})

# getFeedingLevel -----------------------------------------

test_that("getFeedingLevel for MizerParams", {
    fl <- getFeedingLevel(params, n, n_full)
    # test dim
    expect_identical(dim(fl), c(no_sp, no_w))
    expect_identical(dimnames(fl), dimnames(params@initial_n))
    # A crap test - just returns what's already in the function
    encounter <- getEncounter(params, n = n, n_pp = n_full)
    f <- encounter / (encounter + params@intake_max)
    expect_equal(fl, f, ignore_attr = TRUE)
    # test value
    # expect_known_value(fl, "values/getFeedingLevel")
    # expect_snapshot(round(fl, 5)) # round to take into account different rounding errors depending on OS
    expect_snapshot_value(drop_params(fl), style = 'json2', tolerance = 1e-5) # round to take into account different rounding errors depending on OS
})

test_that("getFeedingLevel for MizerSim", {
    time_range <- 1:2
    fl <- getFeedingLevel(sim, time_range = time_range)
    expect_length(dim(fl), 3)
    # because t_save is 0.5, there should be 3 time steps in the range 1:2
    expect_equal(dim(fl), c(3, dim(params@initial_n)))
    expect_identical(dimnames(fl)$sp, dimnames(params@initial_n)$sp)
    expect_identical(dimnames(fl)$w, dimnames(params@initial_n)$w)
    time_range <- 2
    expect_length(dim(getFeedingLevel(sim, time_range = time_range)), 3)
    expect_equal(
        getFeedingLevel(sim, time_range = time_range)[1, , ],
        getFeedingLevel(sim@params, sim@n[as.character(time_range), , ],
                        sim@n_pp[as.character(time_range), ]),
        ignore_attr = TRUE
    )
})

test_that("getFeedingLevel passes correct time", {
    # Here we will check that when getFeedingLevel() is called with
    # a sim object, it passes the correct values of t and n at each time step.
    # To do this we replace mizerFeedingLevel() with a simpler function that
    # just returns t * n
    time_range <- 1:2
    time_elements <- get_time_elements(sim, time_range)
    times <- as.numeric(dimnames(sim@effort)$time[time_elements])
    e <- globalenv() # We need to define the following functions in the
    # global environment so that mizer can find them
    e$testFeedingLevel <- function(params, n, t, ...) {
        n * t
    }
    sim@params <- setRateFunction(sim@params, "FeedingLevel",
                                  "testFeedingLevel")
    expect_equal(getFeedingLevel(sim, time_range = time_range),
                 sweep(sim@n[time_elements, , ], 1, times, "*"),
                 ignore_attr = TRUE)
})

test_that("getFeedingLevel is independent of volume", {
    fl <- getFeedingLevel(params)
    fl_r <- getFeedingLevel(params_r)
    expect_equal(fl, fl_r, ignore_attr = "params")
})

test_that("species-size rate getters work for MizerSim", {
    time_range <- 1:2
    time_elements <- get_time_elements(sim, time_range)
    time_idx <- which(time_elements)[1]
    t <- as.numeric(dimnames(sim@n)$time[[time_idx]])
    n_slice <- array(sim@n[time_idx, , ], dim = dim(sim@n)[2:3])
    dimnames(n_slice) <- dimnames(sim@n)[2:3]
    n_other_slice <- sim@n_other[time_idx, ]
    names(n_other_slice) <- dimnames(sim@n_other)$component
    n_pp_slice <- sim@n_pp[time_idx, ]
    effort_slice <- sim@effort[time_idx, ]

    rates <- list(
        list(
            sim = function() getEncounter(sim, time_range = time_range),
            params = function() getEncounter(sim@params, n = n_slice,
                                             n_pp = n_pp_slice,
                                             n_other = n_other_slice, t = t)
        ),
        list(
            sim = function() getEReproAndGrowth(sim, time_range = time_range),
            params = function() getEReproAndGrowth(sim@params, n = n_slice,
                                                   n_pp = n_pp_slice,
                                                   n_other = n_other_slice,
                                                   t = t)
        ),
        list(
            sim = function() getMort(sim, time_range = time_range),
            params = function() getMort(sim@params, n = n_slice,
                                        n_pp = n_pp_slice,
                                        n_other = n_other_slice,
                                        effort = effort_slice, t = t)
        ),
        list(
            sim = function() getERepro(sim, time_range = time_range),
            params = function() getERepro(sim@params, n = n_slice,
                                          n_pp = n_pp_slice,
                                          n_other = n_other_slice, t = t)
        ),
        list(
            sim = function() getEGrowth(sim, time_range = time_range),
            params = function() getEGrowth(sim@params, n = n_slice,
                                           n_pp = n_pp_slice,
                                           n_other = n_other_slice, t = t)
        ),
        list(
            sim = function() getDiffusion(sim, time_range = time_range),
            params = function() getDiffusion(sim@params, n = n_slice,
                                             n_pp = n_pp_slice,
                                             n_other = n_other_slice, t = t)
        ),
        list(
            sim = function() getFlux(sim, time_range = time_range),
            params = function() getFlux(sim@params, n = n_slice,
                                        n_pp = n_pp_slice,
                                        n_other = n_other_slice, t = t)
        )
    )

    for (rate in rates) {
        sim_rate <- rate$sim()
        expect_s3_class(sim_rate, "ArrayTimeBySpeciesBySize")
        expect_equal(dim(sim_rate), c(sum(time_elements), dim(params@initial_n)))
        expect_equal(sim_rate[1, , ], rate$params(), ignore_attr = TRUE)
    }
})

test_that("getCriticalFeedingLevel matches metab over intake_max times alpha", {
    expected <- params@metab / params@intake_max / params@species_params$alpha
    expect_equal(getCriticalFeedingLevel(params), expected,
                 ignore_attr = c("value_name", "units", "type", "class",
                                 "params", "representation"))
})

# getPredRate -------------------------------------------------------------

test_that("getPredRate for MizerParams", {
    pr <- getPredRate(params, n, n_full)
    expect_s3_class(pr, "ArraySpeciesBySize")
    # test dim
    expect_identical(dim(pr), c(no_sp, no_w_full))
    expect_identical(dimnames(pr)$sp, dimnames(params@initial_n)$sp)
    expect_identical(dimnames(pr)$w_prey, as.character(signif(params@w_full, 3)))
    # test value
    # expect_known_value(pr, "values/getPredRate")
    # expect_snapshot(pr)
    expect_snapshot_value(drop_params(pr), style = 'json2', tolerance = 1e-5) # round to take into account different rounding errors depending on OS
})

test_that("getPredRate for MizerSim", {
    time_range <- 1:2
    time_elements <- get_time_elements(sim, time_range)
    time_idx <- which(time_elements)[1]
    t <- as.numeric(dimnames(sim@n)$time[[time_idx]])
    n_slice <- array(sim@n[time_idx, , ], dim = dim(sim@n)[2:3])
    dimnames(n_slice) <- dimnames(sim@n)[2:3]
    n_other_slice <- sim@n_other[time_idx, ]
    names(n_other_slice) <- dimnames(sim@n_other)$component
    n_pp_slice <- sim@n_pp[time_idx, ]

    pr <- getPredRate(sim, time_range = time_range)
    expect_s3_class(pr, "ArrayTimeBySpeciesBySize")
    expect_equal(dim(pr), c(sum(time_elements), no_sp, no_w_full))
    expect_identical(dimnames(pr)$sp, dimnames(params@initial_n)$sp)
    expect_identical(dimnames(pr)$w_prey,
                     as.character(signif(params@w_full, 3)))
    expect_equal(
        pr[1, , ],
        getPredRate(sim@params, n = n_slice, n_pp = n_pp_slice,
                    n_other = n_other_slice, t = t),
        ignore_attr = TRUE
    )
})

test_that("getPredRate is independent of volume", {
    pr <- getPredRate(params)
    pr_r <- getPredRate(params_r)
    expect_equal(pr, pr_r, ignore_attr = "params")
})

# getPredMort -------------------------------------------------------------------

test_that("getPredMort for MizerParams", {
    # Randomize selectivity and catchability for proper test
    set.seed(0)
    params@catchability[] <-
        runif(prod(dim(params@catchability)), min = 0, max = 1)
    params@selectivity[] <-
        runif(prod(dim(params@selectivity)), min = 0, max = 1)
    m <- getPredMort(params, n, n_full)
    # expect_known_value(m, "values/getPredMort")
    # expect_snapshot(m)
    expect_snapshot_value(drop_params(m), style = 'json2', tolerance = 1e-5) # round to take into account different rounding errors depending on OS

    # Look at numbers in a single prey
    w_offset <- no_w_full - no_w
    m2temp <- rep(NA, no_w)
    pred_rate <- getPredRate(params, n, n_full)
    sp <- runif(1, min = 1, max = no_sp)
    for (i in 1:no_w) {
        m2temp[i] <- sum(params@interaction[, sp] * pred_rate[, w_offset + i])
    }
    expect_equal(m2temp, as.numeric(m[sp, ]))
})

test_that("getPredMort for MizerSim", {
    time_range <- 1:2
    expect_length(dim(getPredMort(sim, time_range = time_range)), 3)
    time_range <- 2
    expect_length(dim(getPredMort(sim, time_range = time_range)), 2)
    ##expect_that(getPredMort(sim, time_range=time_range), equals(getPredMort(sim@params, sim@n[as.character(time_range),,], sim@n_pp[as.character(time_range),])))
    aq1 <- getPredMort(sim, time_range = time_range)
    aq2 <- getPredMort(sim@params, sim@n[as.character(time_range), , ],
                 sim@n_pp[as.character(time_range), ])

    ttot <- 0
    for (i in seq_len(dim(aq1)[1])) {
        ttot <- ttot + sum(aq1[i, ] != aq2[i, ])
    }

    expect_equal(ttot, 0)
})

test_that("getPredMort passes correct time", {
    # Here we will check that when getPredMort() is called with
    # a sim object, it passes the correct values of t and n at each time step.
    # To do this we replace mizerFeedingLevel() with a simpler function that
    # just returns t * n
    times <- as.numeric(dimnames(sim@effort)$time)
    e <- globalenv() # We need to define the following functions in the
    # global environment so that mizer can find them
    e$testPredMort <- function(params, n, n_pp, t, ...) {
        n * t
    }
    sim@params <- setRateFunction(sim@params, "PredMort",
                                  "testPredMort")
    expect_equal(unname(getPredMort(sim)),
                 unname(sweep(sim@n, 1, times, "*")),
                 ignore_attr = TRUE)
})

test_that("interaction is right way round in getPredMort function", {
    sp_name <- NS_species_params_gears_small$species[2]
    inter_small[, sp_name] <- 0  # species not eaten by anything
    params <- newMultispeciesParams(NS_species_params_gears_small, inter_small, info_level = 0)
    m2 <- getPredMort(params, get_initial_n(params), params@cc_pp)
    expect_true(all(m2[sp_name, ] == 0))
})

test_that("getPredMort is independent of volume", {
    pr <- getPredMort(params)
    pr_r <- getPredMort(params_r)
    expect_equal(pr, pr_r, ignore_attr = "params")
})

# getResourceMort ---------------------------------------------------------

test_that("getResourceMort", {
    m2 <- getResourceMort(params, n, n_full)
    # test dim
    expect_length(m2, no_w_full)
    # Check number in final prey size group
    m22 <- colSums(getPredRate(params, n, n_full))
    expect_equal(m22, m2, ignore_attr = TRUE)
    m2b1 <- c(getResourceMort(params, n, n_full))
    # test value
    # expect_known_value(m2b1, "values/getResourceMort")
    # expect_snapshot(m2b1)
    expect_snapshot_value(m2b1, style = 'json2', tolerance = 1e-5) # round to take into account different rounding errors depending on OS
})

test_that("getResourceMort is independent of volume", {
    pm <- getResourceMort(params)
    pm_r <- getResourceMort(params_r)
    expect_equal(pm, pm_r, ignore_attr = TRUE)
})

test_that("getResourceMort and getZ aliases delegate exactly", {
    expect_equal(
        getResourceMort(params, n, n_full),
        colSums(getPredRate(params, n, n_full)),
        ignore_attr = TRUE
    )
    expect_identical(
        getZ(params, n = n, n_pp = n_full, effort = 0.2),
        getMort(params, n = n, n_pp = n_full, effort = 0.2)
    )
})

# getFmortGear ------------------------------------------------------------

test_that("getFmortGear", {
    # Two methods:
    # MizerParams + numeric
    # MizerParams + matrix
    # Randomize selectivity and catchability for proper test
    set.seed(0)
    params@catchability[] <-
        runif(prod(dim(params@catchability)), min = 0, max = 1)
    params@selectivity[] <-
        runif(prod(dim(params@selectivity)), min = 0, max = 1)
    no_gear <- dim(params@catchability)[1]
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    # Single numeric
    effort_num1 <- runif(1, min = 0.1, max = 1)
    # Numeric vector
    effort_num2 <- runif(no_gear, min = 0.1, max = 1)
    # Matrix (or 2D  array) - here with 7 timesteps
    effort_mat <-
        array(runif(no_gear * 7, min = 0.1, max = 1), dim = c(7, no_gear))
    # Call both methods with different effort inputs
    f1 <- getFMortGear(params, effort_num1)
    f2 <- getFMortGear(params, effort_num2)
    f3 <- getFMortGear(params, effort_mat)
    # Check dimnames are right
    expect_named(dimnames(f1), c("gear", "sp", "w"))
    expect_named(dimnames(f2), c("gear", "sp", "w"))
    expect_named(dimnames(f3)[2:4], c("gear", "sp", "w"))
    expect_identical(dim(f3), c(dim(effort_mat)[1], no_gear, no_sp, no_w))
    # check fails if effort is not right size
    bad_effort <- rep(effort_num1, no_gear - 1)
    expect_error(getFMortGear(params, bad_effort))
    # Check contents of output
    widx <- round(runif(1, min = 1, max = no_w))
    sp <- round(runif(1, min = 1, max = no_sp))
    gear <- round(runif(1, min = 1, max = no_gear))
    expect_identical(f1[gear, sp, widx],
                     effort_num1 * params@catchability[gear, sp] *
                         params@selectivity[gear, sp, widx])
    expect_identical(f2[gear, sp, widx],
                     effort_num2[gear] * params@catchability[gear, sp] *
                         params@selectivity[gear, sp, widx])
    expect_identical(f3[, gear, sp, widx],
                     effort_mat[, gear] * params@catchability[gear, sp] *
                         params@selectivity[gear, sp, widx])
    # expect_known_value(f3, "values/getFMortGear")
    # expect_snapshot(f3)
    expect_snapshot_value(f3, style = 'json2', tolerance = 1e-5) # round to take into account different rounding errors depending on OS

    expect_equal(getFMortGear(sim)[1, 1, 1, ],
                 getFMortGear(sim@params, effort = sim@effort[1, ])[1, 1, ])

    time_range <- 1:2
    time_elements <- get_time_elements(sim, time_range)
    f_sim <- getFMortGear(sim, time_range = time_range)
    expected <- plyr::aaply(which(time_elements), 1, function(time_idx) {
        n_slice <- array(sim@n[time_idx, , ], dim = dim(sim@n)[2:3])
        dimnames(n_slice) <- dimnames(sim@n)[2:3]
        n_other_slice <- sim@n_other[time_idx, ]
        names(n_other_slice) <- dimnames(sim@n_other)$component
        getFMortGear(sim@params, effort = sim@effort[time_idx, ],
                     n = n_slice, n_pp = sim@n_pp[time_idx, ],
                     n_other = n_other_slice,
                     t = as.numeric(dimnames(sim@n)$time[[time_idx]]))
    }, .drop = FALSE)
    names(dimnames(expected))[[1]] <- "time"
    expect_equal(f_sim, expected)
})

# getFMort ----------------------------------------------------------------

test_that("getFMort", {
    effort1 <- 0.5
    effort2 <- rep(effort1, no_gear)
    effort3 <- array(effort1, dim = c(7, no_gear))
    f1 <- getFMort(params, effort1)
    f2 <- getFMort(params, effort2)
    f3 <- getFMort(params, effort3)
    # check that length of dims is right
    expect_identical(dim(f1), c(no_sp, no_w))
    expect_identical(dim(f2), c(no_sp, no_w))
    expect_identical(dim(f3), c(dim(effort3)[1], no_sp, no_w))
    # Check dimnames are right
    expect_named(dimnames(f3)[2:3], c("sp", "w"))
    # check fails if effort is not right size
    expect_error(getFMort(params, c(1, 2)))
    # check contents of output
    fmg1 <- getFMortGear(params, effort1)
    fmg2 <- getFMortGear(params, effort2)
    fmg3 <- getFMortGear(params, effort3)
    fmg11 <- array(0, dim = c(no_sp, no_w))
    fmg22 <- array(0, dim = c(no_sp, no_w))
    fmg33 <- array(0, dim = c(dim(effort3)[1], no_sp, no_w))
    for (i in 1:no_gear) {
        fmg11 <- fmg11 + fmg1[i, , ]
        fmg22 <- fmg22 + fmg2[i, , ]
        fmg33 <- fmg33 + fmg3[, i, , ]
    }
    expect_equal(f1, fmg11, ignore_attr = TRUE)
    expect_equal(f2, fmg22, ignore_attr = TRUE)
    expect_equal(f3, fmg33)
    # expect_known_value(f1, "values/getFMort")
    # expect_snapshot(f1)
    expect_snapshot_value(drop_params(f1), style = 'json2', tolerance = 1e-5) # round to take into account different rounding errors depending on OS
})

test_that("getFMort drop argument controls singleton dimensions for MizerSim", {
    single <- project(newMultispeciesParams(NS_species_params_gears_small[3, ], info_level = 0),
                      effort = 1, t_max = 2, progress_bar = FALSE)
    expect_equal(dim(getFMort(single, drop = FALSE)),
                 c(length(getTimes(single)), 1, length(single@params@w)))
    expect_equal(dim(getFMort(single)),
                 c(length(getTimes(single)), length(single@params@w)))
})

test_that("getFMort passes correct time", {
    # Here we will check that when getFMort() calls mizerFMort().
    # it passes the correct time. To implement the test we write simple
    # replacement for mizerFMort() that puts the time
    # argument into the returned array.
    times <- as.numeric(dimnames(sim@effort)$time)
    e <- globalenv() # We need to define the following functions in the
    # global environment so that mizer can find them
    e$testFMort <- function(params, n, t, ...) {
        n * t
    }
    sim@params <- setRateFunction(sim@params, "FMort", "testFMort")
    expect_equal(getFMort(sim),
                 sweep(sim@n, 1, times, "*"),
                 ignore_attr = TRUE)
})

test_that("getFMort passes correct time", {
    # Here we will check that when getFMort() calls getEGrowth().
    # it passes the correct time. To implement the test we write simple
    # replacements for mizerFMort() and mizerEGrowth that put the time
    # argument into the returned array.
    times <- as.numeric(dimnames(sim@effort)$time)
    e <- globalenv() # We need to define the following functions in the
    # global environment so that mizer can find them
    e$testEGrowth <- function(params, n, t, ...) {
        n * t
    }
    e$testFMort <- function(params, e_growth, pred_mort, t, ...) {
        e_growth
    }
    sim@params <- setRateFunction(sim@params, "EGrowth", "testEGrowth")
    sim@params <- setRateFunction(sim@params, "FMort", "testFMort")
    expect_equal(getFMort(sim),
                 sweep(sim@n, 1, times, "*"),
                 ignore_attr = TRUE)

    # Now we do the same for when getFMort() calls getPredMort()
    e$testPredMort <- function(params, n, t, ...) {
        n * t
    }
    e$testFMort <- function(params, e_growth, pred_mort, t, ...) {
        pred_mort
    }
    sim@params <- setRateFunction(sim@params, "PredMort", "testPredMort")
    sim@params <- setRateFunction(sim@params, "FMort", "testFMort")
    expect_equal(getFMort(sim),
                 sweep(sim@n, 1, times, "*"),
                 ignore_attr = TRUE)
})

# getMort --------------------------------------------------------------------

test_that("getMort", {
    no_gear <- dim(params@catchability)[1]
    effort1 <- 0.5
    effort2 <- rep(effort1, no_gear)
    z <- getMort(params, n, n_full, effort = effort2)
    # test dim
    expect_identical(dim(z), c(no_sp, no_w))
    expect_identical(dimnames(z)$sp, dimnames(params@initial_n)$sp)
    expect_identical(dimnames(z)$w, dimnames(params@initial_n)$w)
    # Look at numbers in species 1
    f <- getFMort(params, effort2)
    m2 <- getPredMort(params, n, n_full)
    z1 <- f[1, ] + m2[1, ] + params@species_params$z0[1]
    expect_equal(z1, z[1, ])
    # expect_known_value(z, "values/getMort")
    # expect_snapshot(z)
    expect_snapshot_value(drop_params(z), style = 'json2', tolerance = 1e-5) # round to take into account different rounding errors depending on OS
})

test_that("getMort is independent of volume", {
    m <- getMort(params, effort = 1)
    m_r <- getMort(params_r, effort = 1)
    expect_equal(m, m_r, ignore_attr = "params")
})

test_that("getMort passes state through to getFMort", {
    e <- globalenv()
    e$testFMort <- function(params, n, t, ...) {
        n * t
    }
    e$testMort <- function(params, f_mort, ...) {
        f_mort
    }
    params_f <- setRateFunction(params, "FMort", "testFMort")
    params_f <- setRateFunction(params_f, "Mort", "testMort")
    expect_equal(getMort(params_f, n = n, n_pp = n_full, effort = 1, t = 3),
                 n * 3,
                 ignore_attr = TRUE)
})

test_that("getM2 and getM2Background are aliases", {
    expect_identical(getM2(params, n, n_full), getPredMort(params, n, n_full))
    expect_identical(getM2Background(params, n, n_full),
                     getResourceMort(params, n, n_full))
})

# getEReproAndGrowth ------------------------------------------------------

test_that("getEReproAndGrowth", {
    erg <- getEReproAndGrowth(params, n, n_full)
    # test dim
    expect_identical(dim(erg), c(no_sp, no_w))
    expect_identical(dimnames(erg), dimnames(params@initial_n))
    # Check number in final prey size group
    f <- getFeedingLevel(params, n = n, n_pp = n_full)
    e <-  (f[1, ] * params@intake_max[1, ]) * params@species_params$alpha[1]
    e <- e - params@metab[1, ]
    expect_equal(e, erg[1, ])
    # Can be used with infinite intake_max
    params@intake_max[] <- Inf
    expect_true(!any(is.na(getEReproAndGrowth(params, n = n, n_pp = n_full))))

    erg[erg <= 0] <- 0
    # expect_known_value(erg, "values/getEReproAndGrowth")
    # expect_snapshot(erg)
    expect_snapshot_value(drop_params(erg), style = 'json2', tolerance = 1e-5) # round to take into account different rounding errors depending on OS
})

test_that("getEReproAndGrowth is independent of volume", {
    g <- getEReproAndGrowth(params)
    g_r <- getEReproAndGrowth(params_r)
    expect_equal(g, g_r, ignore_attr = "params")
})

# getERepro ------------------------------------------------------------

test_that("getERepro", {
    es <- getERepro(params, n, n_full)
    # test dim
    expect_identical(dim(es), c(no_sp, no_w))
    expect_identical(dimnames(es), dimnames(params@initial_n))

    e <- getEReproAndGrowth(params, n = n, n_pp = n_full)
    e_repro <- params@psi * e
    e_repro[e_repro <= 0] <- 0
    expect_equal(es, e_repro, ignore_attr = TRUE)
    e_growth <- getEGrowth(params, n, n_full)
    e_growth_man <- e - es
    e_growth_man[e_growth_man <= 0] <- 0
    expect_equal(e_growth, e_growth_man, ignore_attr = TRUE)
    # expect_known_value(es, "values/getERepro")
    # expect_snapshot(es)
    expect_snapshot_value(drop_params(es), style = 'json2', tolerance = 1e-5) # round to take into account different rounding errors depending on OS
})

test_that("getESpawning is an exact alias for getERepro", {
    expect_identical(
        getESpawning(params, n = n, n_pp = n_full),
        getERepro(params, n = n, n_pp = n_full)
    )
})

# getRDI ------------------------------------------------------------------

test_that("getRDI", {
    sex_ratio <- 0.5
    rdi <- getRDI(params, n, n_full)
    # test dim
    expect_length(rdi, no_sp)
    expect_named(rdi, params@species_params$species)
    # test values
    e_repro <- getERepro(params, n = n, n_pp = n_full)
    e_repro_pop <- apply(sweep(e_repro * n, 2, params@dw, "*"), 1, sum)
    rdix <- sex_ratio * (e_repro_pop * params@species_params$erepro) /
        params@w[params@w_min_idx]
    expect_equal(rdix, rdi, tolerance = 1e-15)
    # expect_known_value(rdi, "values/getRDI")
    # expect_snapshot(rdi)
    expect_snapshot_value(rdi, style = 'json2', tolerance = 1e-5) # round to take into account different rounding errors depending on OS
})

test_that("getRDI is proportional to volume", {
    rdi <- getRDI(params)
    rdi_r <- getRDI(params_r)
    expect_equal(rdi * volume, rdi_r)
})

test_that("getRDI for MizerSim", {
    time_range <- 1:2
    time_elements <- get_time_elements(sim, time_range)
    time_idx <- which(time_elements)[1]
    t <- as.numeric(dimnames(sim@n)$time[[time_idx]])
    n_slice <- array(sim@n[time_idx, , ], dim = dim(sim@n)[2:3])
    dimnames(n_slice) <- dimnames(sim@n)[2:3]
    n_other_slice <- sim@n_other[time_idx, ]
    names(n_other_slice) <- dimnames(sim@n_other)$component
    n_pp_slice <- sim@n_pp[time_idx, ]

    rdi <- getRDI(sim, time_range = time_range)
    expect_s3_class(rdi, "ArrayTimeBySpecies")
    expect_equal(dim(rdi), c(sum(time_elements), no_sp))
    expect_identical(dimnames(rdi)$sp, params@species_params$species)
    expect_equal(
        rdi[1, ],
        getRDI(sim@params, n = n_slice, n_pp = n_pp_slice,
               n_other = n_other_slice, t = t),
        ignore_attr = TRUE
    )
})

# getRDD ------------------------------------------------------------------

test_that("getRDD", {
    rdd <- getRDD(params, n, n_full)
    expect_length(rdd, no_sp)
    expect_named(rdd, params@species_params$species)
    rdi <- getRDI(params, n, n_full)
    rdd2 <- getRDD(params, n, n_full, rdi = rdi)
    expect_identical(rdd, rdd2)
    rdd2 <- do.call(params@rates_funcs$RDD,
                    list(rdi = rdi, species_params = params@species_params))
    expect_identical(rdd, rdd2)
    # expect_known_value(rdd, "values/getRDD")
    # expect_snapshot(rdd)
    expect_snapshot_value(rdd, style = 'json2', tolerance = 1e-5) # round to take into account different rounding errors depending on OS
})

test_that("getRDD is proportional to volume", {
    rdd <- getRDD(params)
    rdd_r <- getRDD(params_r)
    expect_equal(rdd * volume, rdd_r)
})

test_that("getRDD for MizerSim", {
    time_range <- 1:2
    time_elements <- get_time_elements(sim, time_range)
    time_idx <- which(time_elements)[1]
    t <- as.numeric(dimnames(sim@n)$time[[time_idx]])
    n_slice <- array(sim@n[time_idx, , ], dim = dim(sim@n)[2:3])
    dimnames(n_slice) <- dimnames(sim@n)[2:3]
    n_other_slice <- sim@n_other[time_idx, ]
    names(n_other_slice) <- dimnames(sim@n_other)$component
    n_pp_slice <- sim@n_pp[time_idx, ]

    rdd <- getRDD(sim, time_range = time_range)
    expect_s3_class(rdd, "ArrayTimeBySpecies")
    expect_equal(dim(rdd), c(sum(time_elements), no_sp))
    expect_identical(dimnames(rdd)$sp, params@species_params$species)
    expect_equal(
        rdd[1, ],
        getRDD(sim@params, n = n_slice, n_pp = n_pp_slice,
               n_other = n_other_slice, t = t),
        ignore_attr = TRUE
    )
})

# getEGrowth --------------------------------------------------------------

test_that("getEGrowth is working", {
    eg1 <- getEGrowth(params, n = n, n_pp = n_full)
    # test dim
    expect_identical(dim(eg1), c(no_sp, no_w))
    expect_identical(dimnames(eg1), dimnames(params@initial_n))
    # expect_known_value(eg1, "values/getEGrowth")
    # expect_snapshot(eg1)
    expect_snapshot_value(drop_params(eg1), style = 'json2', tolerance = 1e-5) # round to take into account different rounding errors depending on OS
})


# getFlux ----

test_that("getFlux works correctly", {
    params <- trait_params_2sp
    # Force different w_min to test zeroing logic
    params@species_params$w_min[2] <- 0.01 
    params@w_min_idx[2] <- which.min(abs(params@w - 0.01))
    
    n <- params@initial_n
    n[] <- 1 # Set n to 1 to make checking easier
    
    t <- 0
    g <- getEGrowth(params, n = n, t = t)
    d <- getDiffusion(params, n = n, t = t)
    dw <- params@dw
    rdd <- getRDD(params, n = n, t = t)
    
    flux <- getFlux(params, n = n, t = t)
    
    # Check dimensions
    expect_equal(dim(flux), dim(n))
    
    # Check zeroing out below w_min_idx
    w_min_idx_2 <- params@w_min_idx[2]
    expect_true(w_min_idx_2 > 1) 
    
    expect_true(all(flux[2, 1:(w_min_idx_2 - 1)] == 0))
    
    # Check boundary condition at w_min_idx
    # flux[i, j_start] = Rdd[i]
    
    # Species 1
    j_start_1 <- params@w_min_idx[1]
    expect_equal(flux[1, j_start_1], rdd[1], ignore_attr = TRUE)
    
    # Species 2
    j_start_2 <- params@w_min_idx[2]
    expect_equal(flux[2, j_start_2], rdd[2], ignore_attr = TRUE)
    
    # Check general calculation for some j > j_start
    # J_{i,j} = g_{i, j-1} N_{i, j-1} - 1/2 * (d_{i, j} N_{i, j} - d_{i, j-1} N_{i, j-1}) / dw_{j-1}
    j <- j_start_2 + 5
    expected_flux_2_j <- g[2, j - 1] * n[2, j - 1] - 0.5 * (d[2, j] * n[2, j] - d[2, j - 1] * n[2, j - 1]) / dw[j - 1]
    
    expect_equal(flux[2, j], expected_flux_2_j, ignore_attr = TRUE)
})

test_that("getFlux uses total diffusion", {
    params <- single_sp_params
    params@use_predation_diffusion <- TRUE
    species <- params@species_params$species[1]

    n <- params@initial_n
    t <- 0
    g <- getEGrowth(params, n = n, t = t)
    d <- getDiffusion(params, n = n, t = t)
    flux <- getFlux(params, n = n, t = t)

    j <- params@w_min_idx[species] + 5
    expected <- g[species, j - 1] * n[species, j - 1] -
        0.5 * (d[species, j] * n[species, j] -
                   d[species, j - 1] * n[species, j - 1]) /
        params@dw[j - 1]

    expect_equal(flux[species, j], expected, ignore_attr = TRUE)
})

test_that("getFlux power argument scales the flux by a power of weight", {
    params <- trait_params_2sp
    n <- params@initial_n

    flux0 <- getFlux(params, n = n)
    flux1 <- getFlux(params, n = n, power = 1)
    flux2 <- getFlux(params, n = n, power = 2)

    # power = 0 is the default
    expect_identical(flux0, getFlux(params, n = n, power = 0))

    expect_equal(flux1[], sweep(flux0[], 2, params@w, "*"),
                 ignore_attr = TRUE)
    expect_equal(flux2[], sweep(flux0[], 2, params@w^2, "*"),
                 ignore_attr = TRUE)

    # units are updated
    expect_equal(attr(flux0, "units"), "1/year")
    expect_equal(attr(flux1, "units"), "g/year")
    expect_equal(attr(flux2, "units"), "g^2/year")

    expect_error(getFlux(params, n = n, power = "a"))
})

test_that("getFlux power argument works for MizerSim", {
    params <- trait_params_2sp
    sim <- project(params, t_max = 1, progress_bar = FALSE)

    flux0 <- getFlux(sim)
    flux1 <- getFlux(sim, power = 1)

    expect_equal(flux1[], sweep(flux0[], 3, params@w, "*"),
                 ignore_attr = TRUE)
    expect_equal(attr(flux1, "units"), "g/year")
})

# getFluxGradient ----

test_that("getFluxGradient returns an ArraySpeciesBySize with the right metadata", {
    fg <- getFluxGradient(NS_params_small)
    expect_s3_class(fg, "ArraySpeciesBySize")
    expect_identical(dim(fg), dim(NS_params_small@initial_n))
    expect_identical(dimnames(fg), dimnames(NS_params_small@metab))
    expect_identical(attr(fg, "value_name"), "Flux gradient")
    expect_identical(attr(fg, "units"), "g^-1/year")
    # The flux gradient is a rate of change of a number density, so it is a
    # density with respect to weight and needs the Jacobian on a length axis
    expect_identical(array_type(fg), "density")
    expect_identical(array_type(getFluxGradient(NS_sim_small)), "density")
    by_weight <- plot(fg, size_axis = "w", return_data = TRUE)
    by_length <- plot(fg, size_axis = "l", return_data = TRUE)
    sp <- NS_params_small@species_params
    sp_idx <- match(by_length$Species, sp$species)
    expect_equal(by_length[[2]],
                 by_weight[[2]] * unname(sp$b[sp_idx]) * by_weight$w /
                     by_length$l)
    expect_identical(plot(fg, size_axis = "l")$scales$get_scales("y")$name,
                     "Flux gradient [cm^-1/year]")
})

test_that("getFluxGradient is the discrete divergence of getFlux in the interior", {
    # (J[j+1] - J[j]) / dw[j] for every bin below the top one. This holds
    # exactly for each flux reconstruction scheme.
    for (scheme in c("upwind", "van_leer", "centred")) {
        params <- NS_params_small
        second_order_w(params) <- c(flux = scheme)
        no_w <- length(params@w)
        no_sp <- nrow(params@species_params)

        fg <- getFluxGradient(params)
        flux <- getFlux(params)
        interior <- (flux[, 2:no_w] - flux[, 1:(no_w - 1)]) /
            matrix(params@dw[1:(no_w - 1)], nrow = no_sp, ncol = no_w - 1,
                   byrow = TRUE)
        expect_equal(unclass(fg)[, 1:(no_w - 1)], interior,
                     ignore_attr = TRUE,
                     info = paste("flux scheme", scheme))
    }
})

test_that("getFluxGradient closes the top bin with a zero phantom bin above", {
    # The upper boundary flux is the scheme's flux formula evaluated with
    # N = 0 in the phantom bin above the grid, so if the density in the top
    # bin is also zero then no flux leaves it and the gradient there is
    # just -J[K] / dw[K].
    for (scheme in c("upwind", "van_leer", "centred")) {
        params <- NS_params_small
        second_order_w(params) <- c(flux = scheme)
        no_w <- length(params@w)

        n <- initialN(params)
        n[, no_w] <- 0
        fg <- getFluxGradient(params, n = n)
        flux <- getFlux(params, n = n)
        expect_equal(unclass(fg)[, no_w], -flux[, no_w] / params@dw[no_w],
                     ignore_attr = TRUE,
                     info = paste("flux scheme", scheme))
    }
})

test_that("getFluxGradient balances mortality at a steady state", {
    # The McKendrick-von Foerster equation is
    #     dN/dt = -dJ/dw - mu N,
    # so at a steady state the flux gradient must cancel the mortality loss.
    params <- NS_params_steady_small
    fg <- getFluxGradient(params)
    loss <- getMort(params) * initialN(params)
    expect_equal(unclass(fg), -unclass(loss), ignore_attr = TRUE,
                 tolerance = 1e-8)
})

test_that("the flux gradient depends on the flux reconstruction scheme", {
    params_up <- NS_params_small
    second_order_w(params_up) <- c(flux = "upwind")
    params_vl <- NS_params_small
    second_order_w(params_vl) <- c(flux = "van_leer")
    expect_false(isTRUE(all.equal(c(getFluxGradient(params_up)),
                                  c(getFluxGradient(params_vl)))))
})

test_that("getFluxGradient for MizerSim matches the MizerParams method", {
    sim <- NS_sim_small
    fg_sim <- getFluxGradient(sim, time_range = c(1, 1))
    expect_s3_class(fg_sim, "ArrayTimeBySpeciesBySize")
    expect_identical(dim(fg_sim), c(1L, dim(sim@params@initial_n)))

    time_idx <- which(get_time_elements(sim, c(1, 1)))[1]
    n_slice <- array(sim@n[time_idx, , ], dim = dim(sim@n)[2:3])
    dimnames(n_slice) <- dimnames(sim@n)[2:3]
    n_other_slice <- sim@n_other[time_idx, ]
    names(n_other_slice) <- dimnames(sim@n_other)$component
    t <- as.numeric(dimnames(sim@n)$time[[time_idx]])

    fg_params <- getFluxGradient(sim@params, n = n_slice,
                                 n_pp = sim@n_pp[time_idx, ],
                                 n_other = n_other_slice, t = t)
    expect_equal(fg_sim[1, , ], fg_params, ignore_attr = TRUE)
})

# Fixed rates let the boundary tests separate transport losses from mortality,
# recruitment feedback and changes in food. Use the ordinary rate dispatch so
# params getters, simulation getters and project() see exactly the same rates.
local_flux_boundary_rates <- function(p, envir = parent.frame()) {
    callbacks <- list(
        flux_boundary_growth = function(params, ...) {
            params@other_params$flux_growth
        },
        flux_boundary_diffusion = function(params, ...) {
            params@other_params$flux_diffusion
        },
        flux_boundary_mort = function(params, ...) params@initial_n * 0,
        flux_boundary_rdd = function(params, ...) {
            rep(0.25, nrow(params@initial_n))
        }
    )
    for (name in names(callbacks)) {
        assign(name, callbacks[[name]], envir = .GlobalEnv)
    }
    withr::defer(rm(list = names(callbacks), envir = .GlobalEnv), envir = envir)
    no_sp <- nrow(p@initial_n)
    p@other_params$flux_growth <- outer(rep(0.1, no_sp), p@w)
    p@other_params$flux_diffusion <- outer(rep(0.1, no_sp), p@w^2)
    p <- setRateFunction(p, "EGrowth", "flux_boundary_growth")
    p <- setRateFunction(p, "Diffusion", "flux_boundary_diffusion")
    p <- setRateFunction(p, "Mort", "flux_boundary_mort")
    setRateFunction(p, "RDD", "flux_boundary_rdd")
}

test_that("flux diagnostics retain upper losses without filling the tail", {
    p0 <- local_flux_boundary_rates(NS_params_small)
    for (scheme in c("upwind", "van_leer", "centred")) {
        for (bin_average in c(FALSE, TRUE)) {
            p <- p0
            second_order_w(p) <- c(flux = scheme, bin_average = bin_average)
            # Different egg sizes also exercise the recruitment boundary.
            p@species_params$w_min[2] <- p@w[3]
            # Let the last species reach the final shared-grid bin, so that
            # the terminal face is covered whatever the fixture's w_max.
            last <- nrow(p@species_params)
            p@species_params$w_max[last] <- max(p@w)
            p <- validParams(p)
            top <- support_top_idx(p)
            cells <- cbind(seq_len(nrow(p@initial_n)), top)
            p@initial_n[] <- 0
            p@initial_n[cells] <- 1 / p@dw[top]
            n <- p@initial_n
            outside <- col(n) > top
            inside <- col(n) >= p@w_min_idx & !outside
            r <- getRates(p)
            raw <- -consumer_residual(p, n, p@initial_n_pp, p@initial_n_other,
                                      p@initial_effort, rates = r)
            outflow <- r$rdd - drop((raw * inside) %*% p@dw)
            expect_true(all(outflow > 0))
            # The last species' outflow must also survive, although getFlux()
            # has no node for that terminal face.
            expect_equal(unname(top[last]), ncol(n))
            flux <- getFlux(p)
            gradient <- getFluxGradient(p)
            sp <- which(top < ncol(n))
            expect_gt(length(sp), 0)
            expect_equal(flux[cbind(sp, top[sp] + 1L)], outflow[sp],
                         ignore_attr = TRUE)
            expect_true(all(gradient[outside] == 0))
            expect_true(all(gradient[col(n) < p@w_min_idx] == 0))
            expect_equal(-gradient[inside], raw[inside], ignore_attr = TRUE)
            expect_equal(-unclass(gradient), steady_rates(p)$dNdt,
                         ignore_attr = TRUE)
            expect_equal(drop(-unclass(gradient) %*% p@dw), r$rdd - outflow,
                         ignore_attr = TRUE)
            eggs <- cbind(seq_len(nrow(n)), p@w_min_idx)
            expect_equal(flux[eggs], r$rdd, ignore_attr = TRUE)
            # An empty retained cell below the top can fill by diffusion.
            expect_true(all(gradient[cbind(seq_len(nrow(n)), top - 1L)] < 0))

            dt <- 1e-6
            sim <- project(p, t_max = dt, dt = dt, t_save = dt,
                           method = "euler", progress_bar = FALSE,
                           info_level = 0)
            observed <- drop((finalN(sim) - n) %*% p@dw) / dt
            expect_equal(observed, r$rdd - outflow, tolerance = 2e-6,
                         ignore_attr = TRUE)
            expect_true(all(finalN(sim)[outside] == 0))
            expect_equal(getFlux(sim, time_range = c(0, 0))[1, , ], flux,
                         ignore_attr = TRUE)
            expect_equal(getFluxGradient(sim, time_range = c(0, 0))[1, , ],
                         gradient, ignore_attr = TRUE)
        }
    }
})

test_that("getFlux excludes tail influx with the solver's original limiter", {
    p0 <- local_flux_boundary_rates(NS_params_small)
    for (scheme in c("upwind", "van_leer", "centred")) {
        for (bin_average in c(FALSE, TRUE)) {
            p <- p0
            second_order_w(p) <- c(flux = scheme, bin_average = bin_average)
            top <- support_top_idx(p)
            n <- p@initial_n * 0
            for (i in seq_len(nrow(n))) {
                n[i, top[i] - 1L] <- 0.25 / p@dw[top[i]]
                n[i, top[i]] <- 1 / p@dw[top[i]]
                if (top[i] < ncol(n)) {
                    n[i, top[i] + 1L] <- 2 / p@dw[top[i]]
                }
            }
            # Preserve the supplied reconstruction history too. Clipping n
            # before evaluating chi would incorrectly change the operator.
            if (scheme != "upwind") {
                attr(n, "chi") <- matrix(0.4, nrow(n), ncol(n))
            }
            p@initial_n <- n
            inside <- col(n) >= p@w_min_idx & col(n) <= top
            r <- getRates(p)
            raw <- -consumer_residual(p, n, p@initial_n_pp, p@initial_n_other,
                                      p@initial_effort, rates = r)
            expected_outflow <- r$rdd - drop((raw * inside) %*% p@dw)
            sp <- which(top < ncol(n))
            boundary <- cbind(sp, top[sp] + 1L)
            flux <- getFlux(p)
            expect_equal(flux[boundary], expected_outflow[sp],
                         ignore_attr = TRUE)
            expect_true(all(flux[boundary] > 0))
            expect_equal(getFlux(p, power = 1)[boundary],
                         flux[boundary] * p@w[top[sp] + 1L],
                         ignore_attr = TRUE)

            # PR #661 retains the algebraic rates of populated excluded bins,
            # but not inflow into empty excluded bins. Fixing the boundary
            # value of getFlux() must not change that chosen convention.
            empty_tail <- col(n) > top & n == 0
            gradient <- getFluxGradient(p)
            expect_equal(-gradient[!empty_tail], raw[!empty_tail],
                         ignore_attr = TRUE)
            expect_true(all(gradient[empty_tail] == 0))
            expect_equal(-unclass(gradient), steady_rates(p)$dNdt,
                         ignore_attr = TRUE)

            # Simulation slices do not store limiter history. Compare getters
            # on exactly the same state without that transient attribute.
            attr(p@initial_n, "chi") <- NULL
            sim <- params_as_sim(p)
            expect_equal(getFlux(sim)[1, , ], getFlux(p), ignore_attr = TRUE)
            expect_equal(getFluxGradient(sim)[1, , ], getFluxGradient(p),
                         ignore_attr = TRUE)
        }
    }
})

test_that("flux diagnostics share the evaluated global Euler boundary policy", {
    p <- local_flux_boundary_rates(NS_params_small)
    second_order_w(p) <- FALSE
    top <- support_top_idx(p)[1]
    p@initial_n[] <- 0
    p@initial_n[1, top] <- 1
    p@other_params$flux_diffusion[] <- 0
    no_diffusion <- getFluxGradient(p)
    boundary_flux <- getFlux(p)[1, top + 1L]
    expect_gt(boundary_flux, 0)
    expect_lt(no_diffusion[1, top + 1L], 0)
    expect_equal(-unclass(no_diffusion), steady_rates(p)$dNdt,
                 ignore_attr = TRUE)
    dt <- 1e-6
    sim <- project(p, t_max = dt, dt = dt, t_save = dt, method = "euler",
                   progress_bar = FALSE, info_level = 0)
    expect_equal(finalN(sim)[1, top + 1L] / dt,
                 -no_diffusion[1, top + 1L], tolerance = 1e-5,
                 ignore_attr = TRUE)

    # Only a different species diffuses, through a custom function rather
    # than ext_diffusion. That still selects the hard boundary for all species.
    expect_true(all(p@ext_diffusion == 0))
    p@other_params$flux_diffusion[2, top] <- 1
    global_diffusion <- getFluxGradient(p)
    expect_equal(unname(global_diffusion[1, top + 1L]), 0)
    expect_equal(global_diffusion[1, top], no_diffusion[1, top])
    expect_equal(getFlux(p)[1, top + 1L], boundary_flux)
    sim <- params_as_sim(p)
    expect_equal(getFluxGradient(sim)[1, , ], global_diffusion,
                 ignore_attr = TRUE)

    for (scheme in c("van_leer", "centred")) {
        second_order_w(p) <- c(flux = scheme)
        top <- support_top_idx(p)[1]
        p@initial_n[] <- 0
        p@initial_n[1, top] <- 1
        p@other_params$flux_diffusion[] <- 0
        expect_gt(getFlux(p)[1, top + 1L], 0)
        expect_equal(unname(getFluxGradient(p)[1, top + 1L]), 0)
    }
})

test_that("getFluxGradient for a simulation uses the Euler boundary policy", {
    # Without diffusion, upwind Euler moves the outflow on into the excluded
    # bin. Predictor-corrector and TR-BDF2 hold that bin at zero instead, as
    # documented, but the gradient does not depend on the simulation's method.
    p <- local_flux_boundary_rates(NS_params_small)
    second_order_w(p) <- FALSE
    p@other_params$flux_diffusion[] <- 0
    top <- support_top_idx(p)[1]
    p@initial_n[] <- 0
    p@initial_n[1, top] <- 1
    n <- p@initial_n
    expected <- getFluxGradient(p)
    expect_lt(expected[1, top + 1L], 0)
    dt <- 1e-6
    for (method in c("euler", "predictor_corrector", "tr_bdf2")) {
        sim <- project(p, t_max = dt, dt = dt, t_save = dt, method = method,
                       progress_bar = FALSE, info_level = 0)
        gradient <- getFluxGradient(sim, time_range = c(0, 0))[1, , ]
        expect_equal(gradient, expected, ignore_attr = TRUE, info = method)
        projected <- (finalN(sim)[1, top + 1L] - n[1, top + 1L]) / dt
        if (method == "euler") {
            expect_equal(unname(-gradient[1, top + 1L]), unname(projected),
                         tolerance = 1e-5)
        } else {
            expect_equal(unname(projected), 0, info = method)
        }
    }
})

test_that("flux diagnostics propagate a non-finite state instead of failing", {
    # Fixed rates keep a non-finite density where it was put. The retained
    # bins never read the excluded neighbour, just as in the solver.
    p0 <- local_flux_boundary_rates(NS_params_small)
    for (scheme in c("upwind", "van_leer")) {
        p <- p0
        second_order_w(p) <- c(flux = scheme)
        top <- support_top_idx(p)
        sp <- which(top < ncol(p@initial_n))
        expect_gt(length(sp), 1)
        retained <- cbind(sp, top[sp])
        tail <- cbind(sp, top[sp] + 1L)
        n <- initialN(p)
        n[tail] <- NaN
        flux <- getFlux(p, n = n)
        gradient <- getFluxGradient(p, n = n)
        expect_true(all(is.finite(flux[retained])))
        expect_true(all(is.finite(flux[tail])))
        expect_true(all(is.finite(gradient[retained])))
        expect_true(all(is.nan(gradient[tail])))
    }

    # With the model's own rates a non-finite density makes every rate
    # non-finite. The diagnostics then return NaN rather than an error, also
    # when the diffusion that selects the boundary policy is itself NaN.
    p <- NS_params_small
    p@use_predation_diffusion <- TRUE
    n <- initialN(p)
    n[] <- NaN
    expect_true(all(is.nan(getFlux(p, n = n))))
    expect_true(all(is.nan(getFluxGradient(p, n = n))))
})
