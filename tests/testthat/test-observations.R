test_that("observation_columns names the columns of each observation type", {
    expect_identical(observation_columns("biomass"),
                     list(to = "biomass",
                          observed = "biomass_observed",
                          cutoff = "biomass_cutoff"))
    expect_identical(observation_columns("number"),
                     list(to = "number",
                          observed = "number_observed",
                          cutoff = "number_cutoff"))
    # Biomass is the default and anything else is rejected
    expect_identical(observation_columns(), observation_columns("biomass"))
    expect_error(observation_columns("yield"))
})

test_that("cutoff_min_w falls back on the smallest weight", {
    params <- NS_params_small
    no_sp <- nrow(params@species_params)

    # No cutoff column at all: the whole size range is counted
    expect_identical(cutoff_min_w(params, "biomass"),
                     rep(min(params@w), no_sp))
    expect_identical(cutoff_min_w(params, "number"),
                     rep(min(params@w), no_sp))

    # A missing entry is filled in, the others are kept
    species_params(params)$biomass_cutoff <- c(10, NA, 20)
    expect_identical(cutoff_min_w(params, "biomass"),
                     c(10, min(params@w), 20))
    # The two types read their own column
    expect_identical(cutoff_min_w(params, "number"),
                     rep(min(params@w), no_sp))
})

test_that("model_observation gives the modelled biomass and number", {
    params <- NS_params_small
    cutoff <- params@w[10]
    species_params(params)$biomass_cutoff <- cutoff
    species_params(params)$number_cutoff <- cutoff

    expect_equal(model_observation(params, "biomass"),
                 getBiomass(params, use_cutoff = TRUE),
                 ignore_attr = TRUE)
    expect_equal(unname(model_observation(params, "number")),
                 unname(getN(params, min_w = cutoff)))

    # Without a cutoff column the whole size range is counted
    no_cutoff <- NS_params_small
    expect_equal(unname(model_observation(no_cutoff, "biomass")),
                 unname(getBiomass(no_cutoff)))
    expect_equal(unname(model_observation(no_cutoff, "number")),
                 unname(getN(no_cutoff)))
})

test_that("model_observation follows the model's quadrature scheme", {
    params <- NS_params_small
    cutoff <- params@w[10]
    species_params(params)$biomass_cutoff <- cutoff
    species_params(params)$number_cutoff <- cutoff
    p2 <- params
    second_order_w(p2) <- c(bin_average = TRUE)

    # The default scheme point-samples the weight at the left bin boundary and
    # cuts the size range at a bin boundary.
    in_range <- params@w >= cutoff
    expect_equal(unname(model_observation(params, "biomass")),
                 unname(rowSums(sweep(params@initial_n, 2,
                                      params@w * params@dw, "*")
                                [, in_range, drop = FALSE])))
    expect_equal(unname(model_observation(params, "number")),
                 unname(rowSums(sweep(params@initial_n, 2, params@dw, "*")
                                [, in_range, drop = FALSE])))

    # With bin averaging on, both the weight and the cutoff mask are averaged
    # over the bin, so both quantities move.
    expect_false(isTRUE(all.equal(unname(model_observation(p2, "biomass")),
                                  unname(model_observation(params, "biomass")))))
    expect_false(isTRUE(all.equal(unname(model_observation(p2, "number")),
                                  unname(model_observation(params, "number")))))
})
