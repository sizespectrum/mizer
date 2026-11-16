test_that("scaleRates with factor = 1 leaves params unchanged", {
    # validParams is called internally, so use it as the base to avoid
    # spurious differences from columns it adds (e.g. D_ext)
    params <- suppressMessages(validParams(NS_params_small))
    expect_identical(scaleRates(params, 1), params)
})

test_that("scaleRates validates factor", {
    params <- NS_params_small
    expect_error(scaleRates(params, 0))
    expect_error(scaleRates(params, -1))
    expect_error(scaleRates(params, "a"))
})

test_that("scaleRates scales all rate slots by factor", {
    params <- NS_params_small
    f <- 3
    scaled <- scaleRates(params, f)

    expect_equal(scaled@search_vol,    params@search_vol    * f)
    expect_equal(scaled@intake_max,    params@intake_max    * f)
    expect_equal(scaled@metab,         params@metab         * f)
    expect_equal(scaled@mu_b,          params@mu_b          * f)
    expect_equal(scaled@ext_encounter, params@ext_encounter * f)
    expect_equal(scaled@ext_diffusion, params@ext_diffusion * f)
    expect_equal(scaled@catchability,  params@catchability  * f)
    expect_equal(scaled@rr_pp,         params@rr_pp         * f)
})

test_that("scaleRates scales species_params columns by factor", {
    params <- NS_params_small
    f <- 3
    scaled <- scaleRates(params, f)
    sp  <- params@species_params
    sp2 <- scaled@species_params

    for (col in c("gamma", "h", "ks", "k", "z0", "E_ext")) {
        if (col %in% names(sp)) {
            expect_equal(sp2[[col]], sp[[col]] * f,
                         label = paste0("species_params$", col))
        }
    }
    # R_max may be Inf by default; Inf * f == Inf
    if ("R_max" %in% names(sp)) {
        expect_equal(sp2$R_max, sp$R_max * f)
    }
})

test_that("scaleRates scales given_species_params columns by factor", {
    params <- NS_params_small
    f <- 3
    scaled <- scaleRates(params, f)
    gsp  <- params@given_species_params
    gsp2 <- scaled@given_species_params

    for (col in c("gamma", "h", "ks", "k", "z0", "R_max")) {
        if (col %in% names(gsp)) {
            expect_equal(gsp2[[col]], gsp[[col]] * f,
                         label = paste0("given_species_params$", col))
        }
    }
})

test_that("scaleRates scales gear_params catchability by factor", {
    params <- NS_params_small
    f <- 3
    scaled <- scaleRates(params, f)
    expect_equal(scaled@gear_params$catchability,
                 params@gear_params$catchability * f)
})

test_that("scaleRates does not add absent species_params columns", {
    params <- NS_params_small
    # z0pre is absent from NS_params_small and should not be invented by scaleRates
    expect_false("z0pre" %in% names(params@species_params))
    scaled <- scaleRates(params, 2)
    expect_false("z0pre" %in% names(scaled@species_params))
})

test_that("scaleRates scales ext_diffusion and D_ext when set", {
    params <- example_params()
    f <- 4
    scaled <- scaleRates(params, f)
    expect_equal(scaled@ext_diffusion, params@ext_diffusion * f)
    if ("D_ext" %in% names(params@species_params)) {
        expect_equal(scaled@species_params$D_ext,
                     params@species_params$D_ext * f)
    }
})

test_that("scaleRates produces time-rescaled dynamics", {
    # Scaling all rates by f is equivalent to running the original model for
    # f times as long with f times larger time steps.
    # project(scaleRates(p, f), t_max=T, dt=dt) should equal
    # project(p, t_max=f*T, dt=f*dt) since each step covers the same
    # biological time.
    f <- 2
    params <- NS_params_small
    params_scaled <- scaleRates(params, f)

    sim_orig   <- project(params,        t_max = 1,   dt = 0.1,  t_save = 1)
    sim_scaled <- project(params_scaled, t_max = 0.5, dt = 0.05, t_save = 0.5)

    # Compare raw arrays to avoid the params attribute carried by ArraySpeciesBySize
    expect_equal(sim_scaled@n[dim(sim_scaled@n)[1], , ],
                 sim_orig@n[dim(sim_orig@n)[1], , ],
                 tolerance = 1e-10)
    expect_equal(sim_scaled@n_pp[dim(sim_scaled@n_pp)[1], ],
                 sim_orig@n_pp[dim(sim_orig@n_pp)[1], ],
                 tolerance = 1e-10)
})

test_that("Rate scaling respects every built-in reproduction law", {
    withr::local_options(mizer_info_level = 0)
    for (average in c(FALSE, TRUE)) {
        for (law in c("BevertonHoltRDD", "RickerRDD", "SheperdRDD",
                      "constantRDD", "noRDD")) {
            p <- NS_params_small
            second_order_w(p) <- c(bin_average = average)
            p <- setReproduction(p, RDD = law)
            if (law == "constantRDD") {
                p <- balanceRDD(p, rdd_target = c(1, 2, 3))
            } else if (law != "noRDD") {
                p <- balanceRDD(p, reproduction_level = 0.4)
            }
            before <- tuning_rate_context(p)
            q <- scaleRates(p, 3)
            after <- tuning_rate_context(q)
            for (field in c("rdi", "rdd", "e_growth", "mort", "diffusion",
                            "resource_mort")) {
                expect_equal(after[[field]], before[[field]] * 3,
                             label = paste(law, average, field))
            }
            expect_equal(reproduction_level(q), reproduction_level(p))
            expect_identical(q$capabilities, p$capabilities)
            expect_identical(q$rates_funcs, p$rates_funcs)
            expect_equal(q$initial_n, p$initial_n)
            expect_equal(q$initial_n_pp, p$initial_n_pp)
            original <- project(p, t_max = 0.3, dt = 0.15, t_save = 0.3)
            scaled <- project(q, t_max = 0.1, dt = 0.05, t_save = 0.1)
            expect_equal(scaled$n[2, , ], original$n[2, , ],
                         tolerance = 1e-10)
            expect_equal(scaled$n_pp[2, ], original$n_pp[2, ],
                         tolerance = 1e-10)
        }
    }
})

test_that("Rate scaling survives rebuilding resource and species defaults", {
    withr::local_options(mizer_info_level = 0)
    for (average in c(FALSE, TRUE)) {
        p <- NS_params_small
        second_order_w(p) <- c(bin_average = average)
        p <- setResource(p, reset = TRUE, balance = FALSE)
        p$given_species_params$z0 <- NULL
        p$given_species_params$gamma <- NULL
        p$given_species_params$ks <- NULL
        p <- setParams(p)
        q <- scaleRates(p, 3)
        expect_equal(q$resource_params$r_pp, p$resource_params$r_pp * 3)
        expect_equal(q$given_species_params$z0, p$species_params$z0 * 3)
        rebuilt <- setResource(q, balance = FALSE, reset = TRUE)
        expect_equal(rebuilt$rr_pp, q$rr_pp)
        rebuilt <- setParams(q)
        for (slot in c("search_vol", "intake_max", "metab", "mu_b")) {
            expect_equal(rebuilt[[slot]], q[[slot]],
                         label = paste(average, slot))
        }
        for (column in c("gamma", "ks", "z0")) {
            expect_equal(rebuilt$species_params[[column]],
                         q$species_params[[column]])
        }
    }
})

test_that("Rate scaling preserves array protection and validates factors", {
    p <- NS_params_small
    slots <- c("search_vol", "intake_max", "metab", "mu_b", "ext_encounter",
               "ext_diffusion", "catchability", "rr_pp")
    for (slot in slots) comment(p[[slot]]) <- paste("manual", slot)
    q <- scaleRates(p, 2)
    for (slot in slots) {
        expect_identical(comment(q[[slot]]), comment(p[[slot]]))
        expect_equal(q[[slot]], p[[slot]] * 2)
    }
    expect_equal(scaleRates(q, 0.5)$species_params, p$species_params)
    for (factor in c(Inf, -Inf, NaN, NA_real_)) {
        expect_error(scaleRates(p, factor))
    }
})

test_that("Rate hooks transform component parameters and retain its state", {
    rlang::local_bindings(
        scaling_component_clock = function(params, before, details, component,
                                           operation, ...) {
            stopifnot(operation == "scale_rates")
            stopifnot(isTRUE(all.equal(params$rr_pp,
                                       before$rr_pp * details$factor)))
            stopifnot(isTRUE(all.equal(params$resource_params$r_pp,
                                       before$resource_params$r_pp *
                                           details$factor)))
            params$other_params[[component]]$turnover <-
                before$other_params[[component]]$turnover * details$factor
            params
        }, .env = globalenv()
    )
    p <- setComponent(NS_params_small, "environment",
        list(stock = 2, temperature = 12, label = "coastal"),
        dynamics_fun = "constant_other",
        component_params = list(turnover = 0.5),
        transformations = list(scale_rates = list(
            action = "update", transform_fun = "scaling_component_clock"
        ))
    )
    q <- scaleRates(p, 4)
    expect_equal(q$other_params$environment$turnover, 2)
    expect_identical(q$initial_n_other, p$initial_n_other)
    expect_identical(q$capabilities, p$capabilities)
    p <- setComponent(p, "unsupported", "intensive",
                      dynamics_fun = "constant_other")
    before <- p
    expect_error(scaleRates(p, 2), "unsupported|support|transformation")
    expect_identical(p, before)
    expect_identical(scaleRates(p, 1), p)
})

test_that("Custom resource clocks require explicit transformation support", {
    rlang::local_bindings(
        scaling_custom_resource = function(params, n_pp, ...) n_pp,
        scaling_resource_clock = function(params, before, details, ...) {
            params$resource_params$period <-
                before$resource_params$period / details$factor
            params
        }, .env = globalenv()
    )
    p <- setResource(NS_params_small,
        resource_dynamics = "scaling_custom_resource", balance = FALSE
    )
    p$resource_params$period <- 12
    before <- p
    expect_error(scaleRates(p, 2), "resource.*support|support.*resource")
    expect_identical(p, before)
    p <- setProcess(p, "resource", transformations = list(
        scale_rates = list(
            action = "update", transform_fun = "scaling_resource_clock",
            paths = list(c("resource_params", "period"))
        )
    ))
    q <- scaleRates(p, 3)
    expect_equal(q$resource_params$period, 4)
    expect_equal(q$rr_pp, p$rr_pp * 3)
    expect_equal(q$resource_params$r_pp, p$resource_params$r_pp * 3)
    expect_identical(q$initial_n_pp, p$initial_n_pp)
})
