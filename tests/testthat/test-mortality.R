test_that("Built-in causes recover the total under both quadratures", {
    for (average in c(FALSE, TRUE)) {
        p <- NS_params_small
        second_order_w(p) <- c(bin_average = average)
        result <- getMortByCause(p)
        expect_s3_class(result, "MizerMortalityBreakdown")
        expect_true(result$complete)
        expect_equal(result$total, getMort(p))
        expect_named(result$causes, c("predation", "fishing", "external"))
        expect_equal(
            as.numeric(result$causes$predation),
            as.numeric(getPredMort(p))
        )
        expect_equal(as.numeric(result$causes$fishing), as.numeric(getFMort(p)))
        expect_equal(as.numeric(result$causes$external), as.numeric(p$mu_b))
        expect_equal(
            as.numeric(Reduce(`+`, result$causes)),
            as.numeric(result$total)
        )
        expect_equal(as.numeric(result$residual), rep(0, length(p$initial_n)),
            tolerance = 1e-12
        )
        expect_equal(nrow(result$unavailable), 0)
        expect_identical(result$metadata$id, names(result$causes))
        expect_output(print(result), "complete")
    }
})

test_that("Component and free mortality causes retain names and styles", {
    rlang::local_bindings(
        mortality_component = function(params, n_other, component, t, ...) {
            params$mu_b * 0 + n_other[[component]] + t
        },
        mortality_free = function(params, multiplier = 1, ...) {
            params$mu_b * 0 + multiplier * 0.1
        }, .env = globalenv()
    )
    p <- setComponent(NS_params_small, "pool", 0.2,
        dynamics_fun = "constant_other",
        mort_fun = "mortality_component", colour = "orange",
        linetype = "dashed"
    )
    other_mort(p)$infection <- "mortality_free"
    result <- getMortByCause(p, t = 2, multiplier = 3)
    expect_true(result$complete)
    expect_named(result$causes, c(
        "predation", "fishing", "external",
        "component:pool", "contribution:infection"
    ))
    expect_equal(
        as.numeric(result$causes[["component:pool"]]),
        rep(2.2, length(p$initial_n))
    )
    expect_equal(
        as.numeric(result$causes[["contribution:infection"]]),
        rep(0.3, length(p$initial_n))
    )
    style <- result$metadata[result$metadata$id == "component:pool", ]
    expect_equal(style$colour, "orange")
    expect_equal(style$linetype, "dashed")
    expect_equal(result$total, getMort(p, t = 2, multiplier = 3))
})

test_that("Constituent overrides retain the supported outer decomposition", {
    rlang::local_bindings(
        mortality_predation = function(params, ...) params$mu_b * 0 + 0.3,
        mortality_fishing = function(params, ...) params$mu_b * 0 + 0.7,
        .env = globalenv()
    )
    p <- setRateFunction(NS_params_small, "PredMort", "mortality_predation")
    p <- setRateFunction(p, "FMort", "mortality_fishing")
    result <- getMortByCause(p)
    expect_true(result$complete)
    expect_equal(
        as.numeric(result$causes$predation),
        rep(0.3, length(p$initial_n))
    )
    expect_equal(
        as.numeric(result$causes$fishing),
        rep(0.7, length(p$initial_n))
    )
})

test_that("Unknown outer providers keep their total without attribution", {
    rlang::local_bindings(
        mortality_outer = function(params, ...) {
            rates <- mizerRates(params, ...)
            rates$mort <- rates$mort + 0.4
            rates
        },
        mortality_total = function(params, ...) params$mu_b * 0 + 0.6,
        .env = globalenv()
    )
    for (rate in c("Rates", "Mort")) {
        callback <- if (rate == "Rates") {
            "mortality_outer"
        } else {
            "mortality_total"
        }
        p <- setRateFunction(NS_params_small, rate, callback)
        result <- getMortByCause(p)
        expect_false(result$complete)
        expect_length(result$causes, 0)
        expected <- tuning_rate_context(p)$mort
        expect_equal(as.numeric(result$total), as.numeric(expected))
        expect_equal(as.numeric(result$residual), as.numeric(expected))
        expect_match(result$unavailable$reason, "declared|current provider")
        expect_output(print(result), "incomplete")
    }
})

test_that("Authored diagnostics share one complete rate context", {
    calls <- 0L
    rlang::local_bindings(
        mortality_counted_rates = function(params, ...) {
            calls <<- calls + 1L
            rates <- mizerRates(params, ...)
            rates$mort <- rates$mort + 0.4
            rates
        },
        mortality_authored = function(p, rates, total, process, t) {
            expect_named(rates, c(
                "encounter", "feeding_level", "e", "e_repro",
                "e_growth", "diffusion", "pred_rate",
                "pred_mort", "f_mort", "mort", "rdi", "rdd",
                "resource_mort"
            ))
            expect_equal(as.numeric(total), as.numeric(rates$mort))
            expect_identical(process, "consumers")
            list(
                causes = list(effective = rates$mort), complete = TRUE,
                metadata = data.frame(
                    id = "effective", label = "All causes",
                    colour = "red"
                )
            )
        }, .env = globalenv()
    )
    p <- setRateFunction(NS_params_small, "Rates", "mortality_counted_rates")
    p <- setProcess(p, "consumers", mortality_fun = "mortality_authored")
    calls <- 0L
    result <- getMortByCause(p, t = 4)
    expect_identical(calls, 1L)
    expect_true(result$complete)
    expect_equal(as.numeric(result$causes$effective), as.numeric(result$total))
    expect_equal(result$metadata$label, "All causes")
    expect_equal(result$metadata$colour, "red")
    # An unrelated rate change leaves this contract current; replacing the
    # outer mortality provider invalidates it until explicitly re-adopted.
    q <- p
    q$rates_funcs$Mort <- "mizer::mizerMort"
    expect_false(getMortByCause(q)$complete)
})

test_that("Mortality method chains need explicit diagnostic adoption", {
    rlang::local_bindings(
        projectMort.mortality_extension = function(params, ...) {
            NextMethod() + 0.2
        },
        mortality_method_breakdown = function(params, rates) {
            list(causes = list(combined = rates$mort), complete = TRUE)
        }, .env = globalenv()
    )
    p <- NS_params_small
    p$extensions <- c(mortality_extension = NA_character_)
    class(p) <- c("mortality_extension", "MizerParams")
    expect_false(builtin_mortality_provider(p))
    result <- getMortByCause(p)
    expect_false(result$complete)
    expect_equal(as.numeric(result$total),
                 as.numeric(getMort(NS_params_small)) + 0.2)
    p <- setProcess(p, "consumers",
                    mortality_fun = "mortality_method_breakdown")
    expect_true(getMortByCause(p)$complete)
    # A qualified direct built-in Mort function bypasses the custom method.
    q <- p
    q$rates_funcs$Mort <- "mizer::mizerMort"
    expect_true(builtin_mortality_provider(q))
    expect_false(getMortByCause(q)$complete) # authored support remains stale
    q <- setProcess(q, "consumers",
                    mortality_fun = "builtin_mortality_breakdown")
    expect_true(getMortByCause(q)$complete)
})

test_that("Incomplete diagnostics preserve residuals and reasons", {
    rlang::local_bindings(
        mortality_partial = function(params, rates) {
            list(
                causes = list(fishing = rates$f_mort), complete = FALSE,
                unavailable = data.frame(
                    id = "other", owner = "consumers",
                    reason = "Other causes not supplied"
                )
            )
        }, .env = globalenv()
    )
    p <- setProcess(NS_params_small, "consumers",
                    mortality_fun = "mortality_partial")
    result <- getMortByCause(p, t = 3)
    expect_false(result$complete)
    expect_equal(result$unavailable$time, 3)
    expect_equal(result$unavailable$reason, "Other causes not supplied")
    expect_equal(
        as.numeric(result$residual),
        as.numeric(result$total - result$causes$fishing)
    )
})

test_that("A projectRates override cannot inherit complete attribution", {
    rlang::local_bindings(
        projectRates.mortality_rates_extension = function(params, ...) {
            rates <- NextMethod()
            rates$mort <- rates$mort * 2
            rates
        }, .env = globalenv()
    )
    p <- NS_params_small
    p$extensions <- c(mortality_rates_extension = NA_character_)
    class(p) <- c("mortality_rates_extension", "MizerParams")
    result <- getMortByCause(p)
    expect_false(result$complete)
    expect_length(result$causes, 0)
    expect_equal(
        as.numeric(result$total),
        as.numeric(getMort(NS_params_small)) * 2
    )
})

test_that("Invalid effective totals and cause names are rejected", {
    rlang::local_bindings(
        mortality_invalid_rates = function(params, ...) {
            rates <- mizerRates(params, ...)
            rates$mort[1, 1] <- Inf
            rates
        },
        mortality_reordered = function(params, rates) {
            list(causes = list(all = rates$mort[3:1, ]), complete = TRUE)
        }, .env = globalenv()
    )
    p <- NS_params_small
    p$rates_funcs$Rates <- "mortality_invalid_rates"
    expect_error(getMortByCause(p), "Effective total.*finite numeric")
    p <- setProcess(NS_params_small, "consumers",
        mortality_fun = "mortality_reordered"
    )
    expect_error(getMortByCause(p), "incompatible species or size names")
})

test_that("Complete mismatches and malformed diagnostics are errors", {
    rlang::local_bindings(
        mortality_bad = function(params, rates) {
            list(causes = list(wrong = rates$mort + 0.1), complete = TRUE)
        }, .env = globalenv()
    )
    p <- setProcess(NS_params_small, "consumers",
                    mortality_fun = "mortality_bad")
    expect_error(getMortByCause(p), "disagrees.*time 0")
    expect_true(getMortByCause(p, atol = 0.2)$complete)
    for (bad in list(
        list(causes = list(), complete = NA),
        list(causes = list(a = 1), complete = TRUE),
        list(
            causes = list(a = NS_params_small$mu_b * NA),
            complete = FALSE
        ),
        list(
            causes = list(a = NS_params_small$mu_b),
            complete = FALSE,
            metadata = data.frame(id = "another")
        )
    )) {
        local({
            value <- bad
            rlang::local_bindings(
                mortality_bad = function(...) value,
                .env = globalenv()
            )
            expect_error(getMortByCause(p), "decomposition|Mortality|metadata")
        })
    }
    rlang::local_bindings(
        mortality_bad = function(...) stop("authored failure"),
        .env = globalenv()
    )
    expect_error(getMortByCause(p), "Mortality decomposition.*authored failure")
    expect_error(getMortByCause(p, rtol = Inf), "finite nonnegative")
    expect_error(getMortByCause(p, atol = -1), "finite nonnegative")
})

test_that("Simulation diagnostics use each saved state, time and effort", {
    rlang::local_bindings(
        mortality_weather = function(params, t, ...) params$mu_b * 0 + t,
        .env = globalenv()
    )
    sim <- NS_sim_small
    sim$effort[2, ] <- sim$effort[2, ] * 2
    other_mort(sim$params)$weather <- "mortality_weather"
    result <- getMortByCause(sim, drop = FALSE)
    expect_s3_class(result$total, "ArrayTimeBySpeciesBySize")
    expect_equal(result$total, getMort(sim, drop = FALSE))
    expect_true(result$complete)
    expect_identical(dim(result$total), dim(sim$n))
    for (index in seq_len(dim(sim$n)[1L])) {
        state <- get_sim_rate_slice(sim, index)
        expected <- do.call(
            getMortByCause,
            c(list(object = sim$params), state)
        )
        expect_equal(
            as.numeric(result$total[index, , ]),
            as.numeric(expected$total)
        )
        expect_equal(as.numeric(result$causes[["contribution:weather"]][
            index, ,
        ]), rep(index - 1, length(sim$params$initial_n)))
    }
    selected <- getMortByCause(sim, t = 2)
    expect_s3_class(selected$total, "ArraySpeciesBySize")
    expect_equal(as.numeric(selected$total), as.numeric(result$total[3, , ]))
})

test_that("Simulation completeness is checked at every time", {
    rlang::local_bindings(
        mortality_cancelled_errors = function(params, rates, t) {
            list(causes = list(all = rates$mort + (t - 1.5)), complete = TRUE)
        }, .env = globalenv()
    )
    sim <- NS_sim_small
    sim$params <- setProcess(sim$params, "consumers",
        mortality_fun = "mortality_cancelled_errors"
    )
    expect_error(getMortByCause(sim), "disagrees.*time 0")
})

test_that("Mortality diagnostics upgrade saved S4 states before defaults", {
    legacy <- readRDS(test_path("assets", "simc.1.0.rds"))
    upgraded <- suppressWarnings(suppressMessages(validSim(legacy)))
    first_time <- as.numeric(dimnames(legacy@n)$time[[1L]])
    result <- suppressWarnings(suppressMessages(
        getMortByCause(legacy, t = first_time)
    ))
    expect_true(result$complete)
    expected <- getMortByCause(upgraded, t = first_time)
    expect_equal(as.numeric(result$total), as.numeric(expected$total))
    result <- suppressWarnings(suppressMessages(getMortByCause(legacy@params)))
    expected <- getMortByCause(upgraded$params)
    expect_true(result$complete)
    expect_equal(as.numeric(result$total), as.numeric(expected$total))
})

test_that("Predator refinement preserves effective total mortality", {
    for (second in c(FALSE, TRUE)) {
        p <- NS_params_small
        p$second_order_w[["bin_average"]] <- second
        aggregate <- getMortByCause(p)
        refined <- getMortByCause(p, predation_by_species = TRUE)
        ids <- paste0("predator:", p$species_params$species)
        expect_true(refined$complete)
        expect_true(all(ids %in% names(refined$causes)))
        expect_false("predation" %in% names(refined$causes))
        expect_equal(Reduce(`+`, refined$causes[ids]),
                     aggregate$causes$predation, ignore_attr = TRUE)
        expect_equal(refined$total, aggregate$total)
        rows <- match(ids, refined$metadata$id)
        expect_identical(refined$metadata$label[rows],
                         as.character(p$species_params$species))
        expect_identical(refined$metadata$colour[rows],
                         unname(getColours(p)[p$species_params$species]))
    }
    sim <- getMortByCause(NS_sim_small, predation_by_species = TRUE,
                          drop = FALSE)
    expect_true(sim$complete)
    ids <- paste0("predator:", NS_params_small$species_params$species)
    expect_true(all(ids %in% names(sim$causes)))
    expect_error(getMortByCause(NS_params_small, predation_by_species = NA),
                 "must be TRUE or FALSE")
})

test_that("Changed predation stays aggregate when refinement disagrees", {
    rlang::local_bindings(
        refined_changed_predation = function(params, ...) {
            2 * mizerPredMort(params, ...)
        }, .env = .GlobalEnv
    )
    p <- setRateFunction(NS_params_small, "PredMort",
                         "refined_changed_predation")
    report <- getMortByCause(p, predation_by_species = TRUE)
    expect_true(report$complete)
    expect_true("predation" %in% names(report$causes))
    expect_false(any(startsWith(names(report$causes), "predator:")))
    expect_equal(report$causes$predation, getPredMort(p), ignore_attr = TRUE)
})

test_that("Author helper shares rates and supports outer assemblers", {
    calls <- 0L
    rlang::local_bindings(
        refined_outer_rates = function(params, n, n_pp, n_other, t,
                                       effort, rates_fns) {
            calls <<- calls + 1L
            rates <- mizerRates(params, n = n, n_pp = n_pp, n_other = n_other,
                                t = t, effort = effort, rates_fns = rates_fns)
            rates$mort <- rates$mort + 0.2
            rates
        },
        refined_author = function(params, ...) {
            result <- mizerMortalityBreakdown(params, ...)
            result$causes$authored <- params$initial_n * 0 + 0.2
            result$metadata <- rbind(result$metadata, data.frame(
                id = "authored", label = "Authored", owner = "consumers",
                colour = "red", linetype = "solid"
            ))
            result
        }, .env = .GlobalEnv
    )
    p <- setRateFunction(NS_params_small, "Rates", "refined_outer_rates")
    p <- setProcess(p, "consumers", mortality_fun = "refined_author")
    calls <- 0L
    report <- getMortByCause(p, predation_by_species = TRUE)
    expect_identical(calls, 1L)
    expect_true(report$complete)
    expect_true("authored" %in% names(report$causes))
    expect_true(any(startsWith(names(report$causes), "predator:")))
    p <- setProcess(p, "consumers", mortality_fun = "mizerMortalityBreakdown")
    expect_error(getMortByCause(p, predation_by_species = TRUE),
                 "Complete mortality decomposition disagrees")
})

test_that("Simulation refinement falls back without reevaluating rates", {
    calls <- 0L
    initial_time <- getTimes(NS_sim_small)[1]
    rlang::local_bindings(
        refined_time_predation = function(params, t, ...) {
            mizerPredMort(params, t = t, ...) * if (t == initial_time) 1 else 2
        },
        refined_counted_rates = function(params, ...) {
            calls <<- calls + 1L
            mizerRates(params, ...)
        }, .env = .GlobalEnv
    )
    p <- setRateFunction(NS_sim_small$params, "PredMort",
                         "refined_time_predation")
    p <- setRateFunction(p, "Rates", "refined_counted_rates")
    p <- setProcess(p, "consumers", mortality_fun = "mizerMortalityBreakdown")
    sim <- NS_sim_small
    sim$params <- p
    calls <- 0L
    report <- getMortByCause(sim, predation_by_species = TRUE, drop = FALSE)
    expect_identical(calls, length(getTimes(sim)))
    expect_true(report$complete)
    expect_true("predation" %in% names(report$causes))
    expect_false(any(startsWith(names(report$causes), "predator:")))
    aggregate <- getMortByCause(sim, drop = FALSE)
    expect_equal(as.numeric(report$causes$predation),
                 as.numeric(aggregate$causes$predation), tolerance = 0)
    expect_null(attr(report, "predation_aggregate"))
    first <- getMortByCause(sim, time_range = initial_time,
                            predation_by_species = TRUE)
    expect_true(any(startsWith(names(first$causes), "predator:")))
})
