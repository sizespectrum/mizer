test_that("component biomass results require explicit valid series", {
    spec <- list(colour = "orange", linetype = "dotted")
    rows <- component_biomass_rows(c(carbon = 2, nitrogen = 1),
                                   "component:detritus", spec)
    expect_equal(rows$id, c("component:detritus/carbon",
                            "component:detritus/nitrogen"))
    expect_equal(rows$colour, rep("orange", 2))
    expect_equal(rows$source, rows$id)
    expect_error(component_biomass_rows(2, "component:x", spec),
                 "series names")
    expect_error(component_biomass_rows(c(x = 1, x = 2), "component:x", spec),
                 "series names")
    expect_error(component_biomass_rows(c(x = Inf), "component:x", spec),
                 "finite and non-negative")
    expect_error(component_biomass_rows(c(x = -1), "component:x", spec),
                 "finite and non-negative")
    expect_error(component_biomass_rows(c(x = NA_real_), "component:x", spec),
                 "finite and non-negative")
    expect_error(component_biomass_rows(numeric(), "component:x", spec),
                 "at least one series")
    missing <- data.frame(series = "carbon", value = NA_real_,
                          available = FALSE, reason = "Unobserved")
    expect_false(component_biomass_rows(missing, "component:x", spec)$available)
    missing$reason <- NA_character_
    expect_error(component_biomass_rows(missing, "component:x", spec),
                 "non-empty reason")
    expect_error(component_biomass_rows(list(carbon = 1), "component:x", spec),
                 "named numeric vector")
})

test_that("biomass series align by identifier across saved times", {
    first <- component_biomass_rows(c(carbon = 2, nitrogen = 1),
                                    "component:detritus", list())
    second <- first[2:1, ]
    second$value <- c(3, 4)
    result <- assemble_biomass(list(first, second), c(0, 1), NS_params_small,
                               simulation = TRUE)
    expect_equal(unname(result[2, ]), c(4, 3))
    second$id[1] <- "component:detritus/other"
    expect_error(assemble_biomass(list(first, second), 0:1, NS_params_small),
                 "identifiers changed")
    second <- first
    second$source[1] <- "another-source"
    expect_error(assemble_biomass(list(first, second), 0:1, NS_params_small),
                 "metadata changed")
})

test_that("biomass selectors reject ambiguous or unavailable input names", {
    p <- NS_params_small
    expect_error(prepare_biomass(p, "not-a-selector", list()),
                 "Unknown biomass selector")
    expect_error(prepare_biomass(p, character(), list()), "non-empty")
    expect_error(prepare_biomass(p, NA_character_, list()), "non-empty")
    expect_error(biomass_requested_filters(list(minimum = 1), FALSE),
                 "must be named")
    expect_error(biomass_requested_filters(list(1), FALSE), "must be named")
    expect_equal(biomass_requested_filters(list(min_w = NULL), FALSE), list())
    expect_equal(biomass_requested_filters(list(min_w = 1), TRUE),
                 list(min_w = 1, use_cutoff = TRUE))
    expect_error(biomass_requested_filters(list(), NA), "TRUE or FALSE")
})

test_that("resource biomass filters have a defined scalar meaning", {
    p <- NS_params_small
    n_pp <- p$initial_n_pp
    p$linecolour["Resource"] <- "purple"
    p$linetype["Resource"] <- "dashed"
    expect_equal(resource_biomass(p, n_pp, list())$colour, "purple")
    expect_equal(resource_biomass(p, n_pp, list())$linetype, "dashed")
    expect_false(resource_biomass(p, n_pp, list(min_w = c(1, 2)))$available)
    expect_false(resource_biomass(p, n_pp, list(use_cutoff = TRUE))$available)
    expect_error(resource_biomass(p, n_pp[-1], list()), "full size grid")
    expect_error(resource_biomass(p, n_pp * NA_real_, list()),
                 "finite non-negative")
    expect_error(resource_biomass(p, rep(1e308, length(n_pp)), list()),
                 "finite and non-negative")
    expect_error(resource_biomass(p, n_pp, list(min_w = NA_real_)),
                 "must not be NA")
    expect_error(resource_biomass(p, n_pp, list(min_l = -1)),
                 "length limits must be non-negative")
    p$resource_params$a <- 0.3
    p$resource_params$b <- 2.5
    length_limit <- 0.1
    weight_limit <- 0.3 * length_limit ^ 2.5
    for (averaged in c(FALSE, TRUE)) {
        p$second_order_w[["bin_average"]] <- averaged
        by_length <- resource_biomass(p, n_pp, list(min_l = length_limit))
        by_weight <- resource_biomass(p, n_pp, list(min_w = weight_limit))
        expect_equal(by_length$value, by_weight$value)
    }
    p$resource_params$b <- NaN
    expect_error(resource_biomass(p, n_pp, list(min_l = 1)),
                 "positive finite weight-length")
})
