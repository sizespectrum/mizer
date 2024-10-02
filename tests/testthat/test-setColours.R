# setColours, getColours ----
test_that("setColours and getColours works", {
    params <- NS_params
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
    expect_false(identical(params@time_modified, NS_params@time_modified))
})

# setLinetypes, getLinetypes ----
test_that("setLinetypes and getLinetypes works", {
    params <- NS_params
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
        "The following are not valid lineypes")
    expect_length(getLinetypes(params), no_types + 2)
    expect_identical(getLinetypes(params)[["test"]], "dotted")
    # Expect updated time_modified
    expect_false(identical(params@time_modified, NS_params@time_modified))
})
