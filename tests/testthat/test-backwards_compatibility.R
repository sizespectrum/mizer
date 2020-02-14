context("Backwards compatibility")

# The known values below were calculated with mizer version 1.0.1

test_that("MizerParams() works as in version 1", {
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    expect_known_value(params@search_vol, "values/set_multispecies_model_search_vol")
    expect_known_value(params@intake_max, "values/set_multispecies_model_intake_max")
})

test_that("set_trait_model() works as in version 1", {
    params <- set_trait_model()
    expect_known_value(params@search_vol, "values/set_trait_model_search_vol")
    expect_known_value(params@intake_max, "values/set_trait_model_intake_max")
})

test_that("set_community_model() works as in version 1", {
    params <- set_community_model()
    expect_known_value(params@search_vol, "values/set_community_model_search_vol")
    expect_known_value(params@intake_max, "values/set_community_model_intake_max")
})
