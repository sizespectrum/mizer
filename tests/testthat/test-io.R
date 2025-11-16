test_that("saveParams/readParams round-trip", {
    params <- NS_params
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    expect_invisible(saveParams(params, tmp))
    params2 <- readParams(tmp)
    expect_s4_class(params2, "MizerParams")
    expect_identical(dim(params2@initial_n), dim(params@initial_n))
    expect_identical(dimnames(params2@initial_n), dimnames(params@initial_n))
    expect_identical(params2@species_params$species, params@species_params$species)
})


