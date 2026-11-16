test_that("Beverton-Holt controls preserve effective current recruitment", {
    p <- NS_params_small
    old <- tuning_rate_context(p)$rdd
    q <- setBevertonHolt(p,
        erepro = 2 * p$species_params$erepro,
        info_level = 0
    )
    expect_equal(q$species_params$erepro, 2 * p$species_params$erepro)
    expect_equal(tuning_rate_context(q)$rdd, old)
    expect_equal(tuning_state(q), tuning_state(p))
    q <- setBevertonHolt(p, R_max = 2 * old, info_level = 0)
    expect_equal(q$species_params$R_max, 2 * old)
    expect_equal(tuning_rate_context(q)$rdd, old)
    expect_equal(unname(reproduction_level(q)), rep(0.5, 3))
    q <- setBevertonHolt(p, reproduction_level = 0.4, info_level = 0)
    expect_equal(tuning_rate_context(q)$rdd, old)
    expect_equal(unname(reproduction_level(q)), rep(0.4, 3))
    expect_false(identical(q$time_modified, p$time_modified))
})
test_that(
    "Beverton-Holt named controls and NAs preserve unselected parameters",
    {
        p <- NS_params_small
        p$given_species_params$R_max <- NULL
        old <- tuning_rate_context(p)$rdd
        q <- setBevertonHolt(p,
            R_max = c(Sprat = 2 * old[[1]], Herring = NA),
            info_level = 0
        )
        expect_equal(q$species_params$R_max[[1]], 2 * old[[1]])
        expect_equal(q$species_params[2:3, ], p$species_params[2:3, ])
        expect_equal(q$given_species_params$R_max[[1]], 2 * old[[1]])
        expect_true(all(is.na(q$given_species_params$R_max[-1])))
        expect_equal(tuning_rate_context(q)$rdd, old)
        q <- setBevertonHolt(p, erepro = c(Sprat = 2), info_level = 0)
        expect_equal(q$species_params$erepro[[1]], 2)
        expect_equal(q$species_params[2:3, ], p$species_params[2:3, ])
        expect_identical(setBevertonHolt(p, erepro = rep(NA, 3)), p)
        expect_identical(setBevertonHolt(p, erepro = c(Sprat = NA)), p)
    }
)
test_that(
    "Beverton-Holt retains erepro by default and rejects infeasible values",
    {
        p <- NS_params_small
        q <- setBevertonHolt(p, info_level = 0)
        expect_equal(q$species_params$erepro, p$species_params$erepro)
        expect_equal(tuning_rate_context(q)$rdd, tuning_rate_context(p)$rdd)
        expect_error(setBevertonHolt(p, erepro = 0), "cannot supply")
        expect_error(setBevertonHolt(p, R_max = 1), "strictly greater")
        expect_error(setBevertonHolt(p, erepro = 1, R_max = 2), "only provide")
        expect_error(setBevertonHolt(p, erepro = "x"), "numeric")
        expect_error(setBevertonHolt(p, R_max = 1:4), "vector of length 3")
        expect_error(setBevertonHolt(p, reproduction_level = 1),
                     "smaller than 1")
        expect_error(setBevertonHolt(p, R_factor = 0.5), "greater than 1")
        expect_error(setBevertonHolt(p, typo = 1), "Unsupported")
    }
)
test_that("Beverton-Holt law switching preserves the previous effective rate", {
    p <- setReproduction(NS_params_small, RDD = "RickerRDD", info_level = 0)
    p <- balanceRDD(p, reproduction_level = 0.4, info_level = 0)
    old <- tuning_rate_context(p)$rdd
    q <- setBevertonHolt(p, reproduction_level = 0.6, info_level = 0)
    expect_identical(q$rates_funcs$RDD, "BevertonHoltRDD")
    expect_equal(tuning_rate_context(q)$rdd, old)
    expect_equal(unname(reproduction_level(q)), rep(0.6, 3))
    expect_error(
        setBevertonHolt(p, reproduction_level = c(Sprat = 0.5)),
        "model-wide"
    )
    withr::local_options(mizer_info_level = 0)
    reproduction_level(p) <- 0.5
    expect_identical(p$rates_funcs$RDD, "RickerRDD")
    expect_equal(tuning_rate_context(p)$rdd, old)
    expect_equal(unname(reproduction_level(p)), rep(0.5, 3))
})
test_that("R_factor and setRmax retain their compatibility semantics", {
    p <- NS_params_small
    q <- setBevertonHolt(p, R_factor = 2, info_level = 0)
    expect_equal(unname(reproduction_level(q)), rep(0.5, 3))
    expect_equal(tuning_rate_context(q)$rdd, tuning_rate_context(p)$rdd)
    a <- setRmax(p, reproduction_level = 0.3, info_level = 0)
    b <- setBevertonHolt(p, reproduction_level = 0.3, info_level = 0)
    expect_equal(a$species_params, b$species_params)
    for (control in list(list(reproduction_level = 0), list(R_factor = Inf))) {
        q <- do.call(setBevertonHolt, c(list(p, info_level = 0), control))
        expect_true(all(is.infinite(q$species_params$R_max)))
        expect_equal(tuning_rate_context(q)$rdd, tuning_rate_context(p)$rdd)
    }
})
test_that(
    "Absent species accept explicit parameters and reject positive levels",
    {
        p <- NS_params_small
        p$initial_n[1, ] <- 0
        q <- setBevertonHolt(p, info_level = 0)
        expect_equal(q$species_params, p$species_params)
        q <- setBevertonHolt(p, erepro = c(Sprat = 0.3), info_level = 0)
        expect_equal(q$species_params$erepro[[1]], 0.3)
        expect_equal(q$given_species_params$erepro[[1]], 0.3)
        expect_equal(q$species_params$R_max, p$species_params$R_max)
        q <- setBevertonHolt(p, R_max = c(Sprat = 1e10), info_level = 0)
        expect_equal(q$species_params$R_max[[1]], 1e10)
        expect_equal(q$species_params$erepro, p$species_params$erepro)
        expect_error(
            setBevertonHolt(p, reproduction_level = c(Sprat = 0.5)),
            "positive reproduction level requires positive recruitment"
        )
        q <- setBevertonHolt(p, reproduction_level = 0, info_level = 0)
        expect_equal(q$species_params$R_max[[1]], p$species_params$R_max[[1]])
        expect_true(all(is.infinite(q$species_params$R_max[-1])))
    }
)
test_that("Information levels govern warnings but not feasibility errors", {
    p <- NS_params_small
    expect_warning(setBevertonHolt(p, erepro = 2), "unrealistic value")
    expect_silent(setBevertonHolt(p, erepro = 2, info_level = 0))
    expect_error(
        setBevertonHolt(p, erepro = 0, info_level = 0),
        "cannot supply"
    )
    withr::local_options(mizer_info_level = 0)
    expect_silent(reproduction_level(p) <- 0.8)
    expect_equal(unname(reproduction_level(p)), rep(0.8, 3))
})
test_that("Law-specific level callbacks respect species names", {
    rlang::local_bindings(
        level_test_reordered = function(params, ...) {
            rev(setNames(
                seq_len(nrow(params$species_params)),
                params$species_params$species
            ))
        },
        level_test_invalid = function(params, ...) {
            c(
                Sprat = 1,
                Sprat = 2, Cod = 3
            )
        },
        .env = globalenv()
    )
    p <- setProcess(NS_params_small, "reproduction",
        level_fun = "level_test_reordered"
    )
    expect_equal(
        reproduction_level(p),
        setNames(1:3, p$species_params$species)
    )
    p <- setProcess(p, "reproduction", level_fun = "level_test_invalid")
    expect_error(reproduction_level(p), "exactly once")
})
test_that("Law switching preserves existing equation arrays", {
    p <- setReproduction(NS_params_small, RDD = "RickerRDD", info_level = 0)
    p$psi <- p$psi * 0.5
    before <- tuning_rate_context(p)
    q <- setBevertonHolt(p, reproduction_level = 0.4, info_level = 0)
    after <- tuning_rate_context(q)
    expect_equal(q$psi, p$psi)
    expect_equal(after$rdd, before$rdd)
    expect_equal(after$e_growth, before$e_growth)
    expect_equal(after$mort, before$mort)
    expect_equal(after$diffusion, before$diffusion)
    expect_equal(tuning_state(q), tuning_state(p))
})
