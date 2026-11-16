test_that("Balancing preserves actual recruitment away from steady state", {
    p <- NS_params_small
    before <- tuning_rate_context(p)
    expect_false(isTRUE(all.equal(before$rdd, getRequiredRDD(p))))
    q <- balanceRDD(p, reproduction_level = 0.4, info_level = 0)
    after <- tuning_rate_context(q)
    expect_equal(after$rdd, before$rdd)
    expect_equal(after$e_growth, before$e_growth)
    expect_equal(after$mort, before$mort)
    expect_equal(after$diffusion, before$diffusion)
    expect_equal(tuning_state(q), tuning_state(p))
    expect_equal(unname(reproduction_level(q)), rep(0.4, 3))
})
test_that("Named targets change only selected species", {
    p <- NS_params_small
    old <- tuning_rate_context(p)$rdd
    q <- balanceRDD(p,
        rdd_target = c(Herring = old[[2]] * 1.2),
        info_level = 0
    )
    expected <- old
    expected[[2]] <- expected[[2]] * 1.2
    expect_equal(tuning_rate_context(q)$rdd, expected)
    expect_equal(q$species_params[c(1, 3), ], p$species_params[c(1, 3), ])
    expect_equal(reproduction_level(q), reproduction_level(p))
    expect_identical(balanceRDD(p, rdd_target = NA), p)
    expect_error(balanceRDD(p, rdd_target = c(unknown = 2)), "unique names")
    expect_error(
        balanceRDD(p, rdd_target = c(Sprat = 1, Sprat = 2)),
        "unique names"
    )
    expect_error(balanceRDD(p, rdd_target = -1), "nonnegative")
    expect_error(balanceRDD(p, rdd_target = Inf), "finite")
    expect_error(balanceRDD(p, rdi = 1), "context cannot")
})
test_that("Explicit preservation and controls are strict", {
    p <- NS_params_small
    target <- tuning_rate_context(p)$rdd
    p$given_species_params$erepro <- NULL
    q <- balanceRDD(p, preserve = "erepro", info_level = 0)
    expect_equal(q$species_params$erepro, p$species_params$erepro)
    expect_null(q$given_species_params$erepro)
    expect_error(balanceRDD(p, erepro = 0), "cannot supply")
    expect_error(balanceRDD(p, R_max = target), "strictly greater")
    expect_error(balanceRDD(p, erepro = 1, R_max = 2), "one explicit")
    expect_error(
        balanceRDD(p, preserve = "erepro", erepro = 1),
        "one explicit"
    )
    expect_error(
        balanceRDD(p, strategy = "set_parameters"),
        "requires an explicit"
    )
    expect_error(balanceRDD(p, imaginary = 1), "Unsupported reproduction")
    q <- balanceRDD(p, reproduction_level = 0, info_level = 0)
    expect_true(all(is.infinite(q$species_params$R_max)))
    expect_equal(tuning_rate_context(q)$rdd, target)
})
test_that(
    "Ricker inversion selects ascending and descending limbs explicitly",
    {
        p <- setReproduction(NS_params_small, RDD = "RickerRDD", info_level = 0)
        old <- tuning_rate_context(p)$rdd
        ascending <- balanceRDD(p, reproduction_level = 0.4, info_level = 0)
        descending <- balanceRDD(p,
            reproduction_level = 0.4,
            limb = "descending", info_level = 0
        )
        a <- tuning_rate_context(ascending)
        d <- tuning_rate_context(descending)
        expect_equal(a$rdd, old)
        expect_equal(d$rdd, old)
        expect_true(all(ascending$species_params$ricker_b * a$rdi < 1))
        expect_true(all(descending$species_params$ricker_b * d$rdi > 1))
        expect_equal(reproduction_level(ascending),
                     reproduction_level(descending))
        again <- balanceRDD(descending, info_level = 0)
        expect_equal(again$species_params$erepro,
                     ascending$species_params$erepro)
        peak <- balanceRDD(p, reproduction_level = 1, info_level = 0)
        expect_equal(unname(reproduction_level(peak)), rep(1, 3))
        expect_error(balanceRDD(p, reproduction_level = 1.1), "must not exceed")
        expect_error(
            balanceRDD(p, ricker_b = 2 / (exp(1) * old)),
            "exceeds the maximum"
        )
        expect_error(
            balanceRDD(p, reproduction_level = 0, limb = "descending"),
            "no descending limb"
        )
    }
)
test_that("Shepherd has law-specific levels and finite inverse branches", {
    p <- setReproduction(NS_params_small, RDD = "SheperdRDD", info_level = 0)
    target <- tuning_rate_context(p)$rdd
    a <- balanceRDD(p,
        reproduction_level = 0.4, sheperd_c = 2,
        info_level = 0
    )
    d <- balanceRDD(p,
        reproduction_level = 0.4, sheperd_c = 2,
        limb = "descending", info_level = 0
    )
    ar <- tuning_rate_context(a)
    dr <- tuning_rate_context(d)
    expect_equal(ar$rdd, target)
    expect_equal(dr$rdd, target)
    expect_true(all(a$species_params$sheperd_b * ar$rdi < 1))
    expect_true(all(d$species_params$sheperd_b * dr$rdi > 1))
    expect_equal(unname(reproduction_level(a)), rep(0.4, 3))
    monotone <- balanceRDD(p,
        reproduction_level = 0.4, sheperd_c = 1,
        info_level = 0
    )
    expect_equal(unname(reproduction_level(monotone)), rep(0.4, 3))
    expect_error(balanceRDD(monotone,
                     preserve = "sheperd_b",
                     limb = "descending"
                 ), "no descending limb")
    p$species_params$sheperd_c[] <- 0.5
    p <- setReproduction(p, info_level = 0)
    expect_equal(unname(reproduction_level(p)), rep(0, 3))
    q <- balanceRDD(p, preserve = "sheperd_b", info_level = 0)
    expect_equal(tuning_rate_context(q)$rdd, tuning_rate_context(p)$rdd)
    expect_error(balanceRDD(p, reproduction_level = 0.4), "no informative")
})
test_that("Linear and constant recruitment have their own level semantics", {
    p <- setReproduction(NS_params_small, RDD = "noRDD", info_level = 0)
    target <- tuning_rate_context(p)$rdd
    q <- balanceRDD(p, rdd_target = target * 0.7, info_level = 0)
    expect_equal(tuning_rate_context(q)$rdd, target * 0.7)
    expect_equal(q$species_params$erepro, p$species_params$erepro * 0.7)
    expect_equal(unname(reproduction_level(q)), rep(0, 3))
    expect_error(balanceRDD(p, reproduction_level = 0.2), "noRDD has level 0")
    p <- setReproduction(p, RDD = "constantRDD", info_level = 0)
    q <- balanceRDD(p, rdd_target = c(1, 0, 3), info_level = 0)
    expect_equal(unname(tuning_rate_context(q)$rdd), c(1, 0, 3))
    expect_equal(unname(reproduction_level(q)), c(1, NA, 1))
    expect_error(balanceRDD(q, reproduction_level = 0.5), "level 1")
})
test_that("Zero efficiency is handled from the declared RDI formula", {
    for (average in c(FALSE, TRUE)) {
        p <- NS_params_small
        p$second_order_w$bin_average <- average
        p$species_params$erepro[] <- 0
        p <- setReproduction(p, info_level = 0)
        q <- balanceRDD(p,
            rdd_target = c(Sprat = 1),
            reproduction_level = c(Sprat = 0.5), info_level = 0
        )
        expect_equal(tuning_rate_context(q)$rdd[[1]], 1)
        expect_gt(q$species_params$erepro[[1]], 0)
        expect_equal(unname(q$species_params$erepro[-1]), c(0, 0))
        p$initial_n[1, ] <- 0
        expect_error(
            balanceRDD(p,
                rdd_target = c(Sprat = 1),
                reproduction_level = c(Sprat = 0.5)
            ),
            "no reproductive energy"
        )
    }
})
test_that(
    "Undeclared effective reproduction providers cannot inherit support",
    {
        rlang::local_bindings(
            balance_test_rdi = function(params, ...) getRequiredRDD(params),
            balance_test_rates = function(params, ...) mizerRates(params, ...),
            .env = globalenv()
        )
        p <- setRateFunction(NS_params_small, "RDI", "balance_test_rdi")
        expect_error(balanceRDD(p), "stationary|support")
        expect_error(reproduction_level(p), "reproduction-level support")
        p <- setRateFunction(NS_params_small, "Rates", "balance_test_rates")
        expect_error(balanceRDD(p), "stationary|support")
    }
)
test_that("S3 balancing chains share the verified transaction in tuning", {
    visits <- character()
    rlang::local_bindings(
        balance_test_immigration = function(rdi, species_params, ...) {
            BevertonHoltRDD(rdi, species_params) + 1
        },
        balanceRDD.balance_outer = function(params, ...) {
            visits <<- c(visits, "outer")
            NextMethod()
        },
        balanceRDD.balance_inner = function(params, rdd_target, ...) {
            visits <<- c(visits, "inner")
            NextMethod(rdd_target = rdd_target - 1)
        },
        .env = globalenv()
    )
    p <- NS_params_small
    tuning <- p$capabilities$processes$reproduction$tuning
    tuning$default <- "preserve_erepro"
    p$extensions <- c(
        balance_outer = NA_character_,
        balance_inner = NA_character_
    )
    class(p) <- c("balance_outer", "balance_inner", class(p))
    p <- setRateFunction(p, "RDD", "balance_test_immigration")
    p <- setProcess(p, "consumers")
    p <- setProcess(p, "resource")
    p <- setProcess(p, "reproduction",
        stationary = TRUE,
        equation_rates_fun = "builtin_reproduction_equation_rates",
        tuning = tuning
    )
    before <- tuning_rate_context(p)$rdd
    q <- balanceRDD(p, info_level = 0)
    expect_identical(visits, c("outer", "inner"))
    expect_equal(tuning_rate_context(q)$rdd, before)
    visits <- character()
    plan <- prepare_tuning(p, owners = "reproduction", solver = "project")
    targets <- capture_tuning_targets(p, plan)
    q <- finish_tuning(p, plan, targets)
    expect_identical(visits, c("outer", "inner"))
    expect_equal(tuning_rate_context(q)$rdd, before)
    rlang::local_bindings(
        balanceRDD.balance_inner = function(params, ...) {
            result <- NextMethod()
            result$initial_n[1, 1] <- result$initial_n[1, 1] * 2
            result
        }, .env = globalenv()
    )
    # Re-register the changed method chain before exercising its candidate.
    p <- setProcess(p, "consumers")
    p <- setProcess(p, "resource")
    p <- setProcess(p, "reproduction")
    expect_error(balanceRDD(p, info_level = 0), "target state")
})
test_that("Declared coupled reproduction rejects partial transactions", {
    p <- NS_params_small
    tuning <- p$capabilities$processes$reproduction$tuning
    tuning$strategies[[tuning$default]]$partial_species <- FALSE
    p <- setProcess(p, "reproduction", tuning = tuning)
    expect_error(
        require_reproduction_balance_support(p, partial = TRUE),
        "partial species"
    )
    expect_error(balanceRDD(p, rdd_target = c(Sprat = 2)), "partial species")
})
test_that("Scalar inverse formulas handle extreme and boundary values", {
    for (limb in c("ascending", "descending")) {
        x <- ricker_inverse(1e-100, 1e-100, limb, "test")
        expect_equal(log(x) - 1e-100 * x, log(1e-100), tolerance = 1e-10)
        x <- shepherd_inverse(0.1, 1, 2, limb, "test")
        expect_equal(x / (1 + x^2), 0.1)
    }
    expect_equal(ricker_inverse(1 / exp(1), 1, "ascending", "test"), 1)
    expect_equal(shepherd_inverse(0.5, 1, 2, "descending", "test"), 1)
    expect_equal(shepherd_inverse(0.2, 1, 1, "ascending", "test"), 0.25)
    expect_error(ricker_inverse(0, 1, "descending", "test"), "finite")
    expect_error(shepherd_inverse(0, 1, 2, "descending", "test"), "finite")
    expect_error(shepherd_inverse(1, 1, 2, "ascending", "test"), "maximum")
})
test_that("A zero target is attainable on an unconstrained linear limit", {
    for (law in c("BevertonHoltRDD", "RickerRDD", "SheperdRDD")) {
        p <- setReproduction(NS_params_small, RDD = law, info_level = 0)
        p <- balanceRDD(p, reproduction_level = 0, info_level = 0)
        q <- balanceRDD(p, rdd_target = c(Sprat = 0), info_level = 0)
        expect_equal(tuning_rate_context(q)$rdd[[1]], 0)
        expect_equal(q$species_params$erepro[[1]], 0)
        expect_equal(q$species_params[-1, ], p$species_params[-1, ])
    }
})

test_that("Qualified built-in reproduction retains declared balancing", {
    p <- setRateFunction(NS_params_small, "RDD", "mizer::BevertonHoltRDD")
    expect_identical(builtin_reproduction_law(p), "BevertonHoltRDD")
    q <- balanceRDD(p, reproduction_level = 0.4, info_level = 0)
    expect_equal(tuning_rate_context(q)$rdd, tuning_rate_context(p)$rdd)
    expect_equal(unname(reproduction_level(q)), rep(0.4, 3))
})


test_that("Built-in reproduction definitions are reusable serialisable data", {
    for (law in c("BevertonHoltRDD", "RickerRDD", "SheperdRDD",
                  "constantRDD", "noRDD")) {
        p <- setReproduction(NS_params_small, RDD = law, info_level = 0)
        definition <- mizerReproductionDefinition(p)
        expect_named(definition, c("stationary", "equation_rates_fun",
                                    "tuning", "level_fun"))
        expect_true(tuning_data_valid(definition))
        expect_equal(definition$tuning, builtin_reproduction_tuning(p))
        q <- setProcess(p, "reproduction", definition = definition)
        expect_true(tuning_binding_current(
            q, q$capabilities$processes$reproduction, "reproduction"))
        expect_true(level_binding_current(
            q, q$capabilities$processes$reproduction))
    }
})

test_that("Definition reuse does not automatically endorse custom routes", {
    p <- NS_params_small
    p$rates_funcs$Mort <- "author_mortality"
    before <- p
    definition <- mizerReproductionDefinition(p)
    expect_identical(p, before)
    expect_null(default_provider_tuning(p, "reproduction")$tuning)
    expect_type(definition$tuning, "list")
    expect_false(tuning_binding_current(
        p, p$capabilities$processes$reproduction, "reproduction"))
    p$rates_funcs$RDD <- "mizer::BevertonHoltRDD"
    expect_identical(mizerReproductionDefinition(p)$level_fun,
                     "reproduction_level_BevertonHoltRDD")
    p$rates_funcs$RDD <- "author_recruitment"
    expect_error(mizerReproductionDefinition(p), "No built-in reproduction")
})
