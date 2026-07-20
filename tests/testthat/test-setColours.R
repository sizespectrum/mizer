# setColours, getColours ----
test_that("setColours and getColours works", {
    params <- NS_params_small
    no_col <- length(getColours(params))
    # nothing changes when setting the same colours
    expect_identical(setColours(params, getColours(params)), params)
    # set new entry
    params <- setColours(params, list("test" = "orange"))
    expect_length(getColours(params), no_col + 1)
    expect_identical(getColours(params)[["test"]], "orange")
    # overwrite existing and set new
    params <- setColours(params, list(test = "blue", test2 = "orange"))
    expect_length(getColours(params), no_col + 2)
    expect_identical(getColours(params)[["test"]], "blue")
    expect_identical(getColours(params)[["test2"]], "orange")
    # NAs are ignored
    params <- setColours(params, list(test = NA))
    expect_identical(getColours(params)[["test"]], "blue")
    # Invalid colours are ignored and trigger warning
    expect_warning(
        params <- setColours(params, list(test = "igit", test3 = "igitigit")),
        "The following are not valid colour values and will be ignored: igit, igitigit")
    expect_length(getColours(params), no_col + 2)
    expect_identical(getColours(params)[["test"]], "blue")
    # Expect updated time_modified
    expect_false(identical(params@time_modified, NS_params_small@time_modified))
})

test_that("setColours updates species_params and given_species_params", {
    params <- setColours(NS_params_small, list(Cod = "black"))
    expect_identical(species_params(params)["Cod", "linecolour"], "black")
    expect_identical(given_species_params(params)["Cod", "linecolour"], "black")
    expect_identical(params@linecolour[["Cod"]], "black")
    # Non-species names are not added to species_params
    expect_null(species_params(params)[["test"]])
    params <- setColours(params, list(test = "orange"))
    expect_null(species_params(params)[["test"]])
    # The colour survives a recalculation of the model
    params <- setBevertonHolt(params, erepro = 0.5)
    expect_identical(species_params(params)["Cod", "linecolour"], "black")
})

# setLinetypes, getLinetypes ----
test_that("setLinetypes and getLinetypes works", {
    params <- NS_params_small
    no_types <- length(getLinetypes(params))
    # nothing changes when using existing types
    expect_identical(setLinetypes(params, getLinetypes(params)), params)
    # set new entry
    params <- setLinetypes(params, list("test" = "dashed"))
    expect_equal(length(getLinetypes(params)), no_types + 1)
    expect_identical(getLinetypes(params)[["test"]], "dashed")
    # overwrite existing and set new
    params <- setLinetypes(params, list("test" = "dotted", test2 = "dashed"))
    expect_equal(length(getLinetypes(params)), no_types + 2)
    expect_identical(getLinetypes(params)[["test"]], "dotted")
    expect_identical(getLinetypes(params)[["test2"]], "dashed")
    # NAs are ignored
    params <- setLinetypes(params, list(test = NA))
    expect_identical(getLinetypes(params)[["test"]], "dotted")
    # Invalid linetypes are ignored and trigger warning
    expect_warning(
        params <- setLinetypes(params, list(test = "igit", test3 = "igitigit")),
        "The following are not valid linetypes")
    expect_length(getLinetypes(params), no_types + 2)
    expect_identical(getLinetypes(params)[["test"]], "dotted")
    # Expect updated time_modified
    expect_false(identical(params@time_modified, NS_params_small@time_modified))
})

test_that("setLinetypes updates species_params and given_species_params", {
    params <- setLinetypes(NS_params_small, list(Cod = "dashed"))
    expect_identical(species_params(params)["Cod", "linetype"], "dashed")
    expect_identical(given_species_params(params)["Cod", "linetype"], "dashed")
    expect_identical(params@linetype[["Cod"]], "dashed")
    # Non-species names are not added to species_params
    params <- setLinetypes(params, list(test = "dotted"))
    expect_null(species_params(params)[["test"]])
    # The linetype survives a recalculation of the model
    params <- setBevertonHolt(params, erepro = 0.5)
    expect_identical(species_params(params)["Cod", "linetype"], "dashed")
})
