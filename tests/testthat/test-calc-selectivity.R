test_that("calc_selectivity builds array and applies selectivity functions", {
    spg <- NS_species_params_gears
    spg$sel_func <- "knife_edge"
    spg$knife_edge_size <- 1000
    spg$knife_edge_size[spg$gear == "Industrial"] <- 500
    spg <- spg[, !(names(spg) %in% c("l25", "l50", "a", "b"))]
    params <- newMultispeciesParams(spg, inter, info_level = 0)
    sel <- calc_selectivity(params)
    expect_equal(dim(sel), dim(params@selectivity))
    industrial_species <- spg$species[spg$gear == "Industrial"]
    expect_true(all(sel["Industrial", industrial_species, params@w < 500] == 0))
    expect_true(all(sel["Industrial", industrial_species, params@w >= 500] == 1))
})

test_that("calc_selectivity leaves unspecified gear-species combinations at zero", {
    spg <- NS_species_params_gears
    spg$sel_func <- "knife_edge"
    spg$knife_edge_size <- 1000
    params <- newMultispeciesParams(spg, inter, info_level = 0)

    sel <- calc_selectivity(params)
    missing_species <- setdiff(params@species_params$species,
                               spg$species[spg$gear == "Industrial"])

    expect_true(all(sel["Industrial", missing_species, ] == 0))
    expect_identical(dimnames(sel), dimnames(params@selectivity))
})

test_that("calc_selectivity errors for missing or NA selectivity parameters", {
    spg <- NS_species_params_gears
    spg$sel_func <- "knife_edge"
    spg$knife_edge_size <- 1000
    spg <- spg[, !(names(spg) %in% c("l25", "l50", "a", "b"))]

    params_missing <- newMultispeciesParams(spg, inter, info_level = 0)
    params_missing@gear_params$knife_edge_size <- NULL
    expect_error(calc_selectivity(params_missing),
                 "missing in the gear_params dataframe")

    params_na <- newMultispeciesParams(spg, inter, info_level = 0)
    params_na@gear_params$knife_edge_size[1] <- NA
    expect_error(calc_selectivity(params_na), "Some selectivity parameters are NA")
})
