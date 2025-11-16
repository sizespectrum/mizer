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


