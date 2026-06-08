knife_edge_selectivity_species_params <- NS_species_params_gears_small
knife_edge_selectivity_species_params$sel_func <- "knife_edge"
knife_edge_selectivity_species_params$knife_edge_size <- 1000
knife_edge_selectivity_params <- newMultispeciesParams(
    knife_edge_selectivity_species_params, inter_small, info_level = 0)

knife_edge_selectivity_species_params_no_length <-
    knife_edge_selectivity_species_params[
        , !(names(knife_edge_selectivity_species_params) %in%
                c("l25", "l50", "a", "b"))
    ]
knife_edge_selectivity_params_no_length <- newMultispeciesParams(
    knife_edge_selectivity_species_params_no_length, inter_small, info_level = 0)

test_that("calc_selectivity builds array and applies selectivity functions", {
    spg <- knife_edge_selectivity_species_params_no_length
    spg$knife_edge_size[spg$gear == "Industrial"] <- 500
    params <- newMultispeciesParams(spg, inter_small, info_level = 0)
    sel <- calc_selectivity(params)
    expect_equal(dim(sel), dim(params@selectivity))
    industrial_species <- spg$species[spg$gear == "Industrial"]
    expect_true(all(sel["Industrial", industrial_species, params@w < 500] == 0))
    expect_true(all(sel["Industrial", industrial_species, params@w >= 500] == 1))
})

test_that("calc_selectivity leaves unspecified gear-species combinations at zero", {
    spg <- knife_edge_selectivity_species_params
    params <- knife_edge_selectivity_params

    sel <- calc_selectivity(params)
    missing_species <- setdiff(params@species_params$species,
                               spg$species[spg$gear == "Industrial"])

    expect_true(all(sel["Industrial", missing_species, ] == 0))
    expect_identical(dimnames(sel), dimnames(params@selectivity))
})

test_that("calc_selectivity errors for missing or NA selectivity parameters", {
    params_missing <- knife_edge_selectivity_params_no_length

    params_missing@gear_params$knife_edge_size <- NULL
    expect_error(calc_selectivity(params_missing),
                 "missing in the gear_params dataframe")

    params_na <- knife_edge_selectivity_params_no_length
    params_na@gear_params$knife_edge_size[1] <- NA
    expect_error(calc_selectivity(params_na), "Some selectivity parameters are NA")
})
