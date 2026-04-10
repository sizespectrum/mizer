test_that("setExtEncounter works", {
    params <- NS_params

    # Without ext_encounter argument, recalculates from E_ext (default 0) and n
    p2 <- setExtEncounter(params)
    zero_encounter <- p2@mu_b # to get right dimnames
    zero_encounter[] <- 0
    expect_identical(zero_encounter, p2@ext_encounter)

    # supplying ext_encounter
    p2 <- setExtEncounter(params, 3 * params@mu_b)
    expect_equal(p2@ext_encounter, 3 * params@mu_b, ignore_attr = TRUE)

    # only ext_encounter changed
    p2@ext_encounter <- params@ext_encounter
    expect_unchanged(p2, params)

    # has updated time_modified
    expect_false(identical(params@time_modified, p2@time_modified))
})

test_that("setExtEncounter uses E_ext species param", {
    params <- NS_params
    species_params(params)$E_ext <- 0.1
    p2 <- setExtEncounter(params, reset = TRUE)
    expected <- sweep(outer(species_params(params)[["n"]],
                            w(params), function(x, y) y^x),
                      1, species_params(params)[["E_ext"]], "*")
    expect_equal(p2@ext_encounter, expected, ignore_attr = TRUE)
})

test_that("reset works on ext_encounter", {
    params <- NS_params
    # Set a custom ext_encounter with a comment
    custom <- params@ext_encounter
    custom[] <- 1
    comment(custom) <- "custom"
    params <- setExtEncounter(params, ext_encounter = custom)
    expect_identical(comment(params@ext_encounter), "custom")

    # reset = TRUE ignores comment and recalculates from species params (E_ext=0)
    p2 <- setExtEncounter(params, reset = TRUE)
    expect_null(comment(p2@ext_encounter))
    expect_equal(p2@ext_encounter, params@ext_encounter * 0, ignore_attr = TRUE)
})

test_that("Comment works on ext_encounter", {
    params <- NS_params
    ext_encounter <- params@ext_encounter
    # comment is stored
    comment(ext_encounter) <- "test"
    params <- setExtEncounter(params, ext_encounter = ext_encounter)
    expect_identical(comment(params@ext_encounter), "test")
    
    # if no comment, previous comment is kept
    comment(ext_encounter) <- NULL
    params <- setExtEncounter(params, ext_encounter = ext_encounter)
    expect_identical(comment(params@ext_encounter), "test")
})

test_that("getExtEncounter works", {
    expect_true(is.ArraySpeciesBySize(getExtEncounter(NS_params)))
    expect_equal(getExtEncounter(NS_params), NS_params@ext_encounter,
                 ignore_attr = TRUE)
})

test_that("setExtEncounter validates dimensions", {
    expect_error(setExtEncounter(NS_params, array(0, dim = c(1, 1))))
})

test_that("Can get and set slot", {
    params <- NS_params
    ext_encounter <- getExtEncounter(params)
    expect_identical(ext_encounter(params), ext_encounter)
    new <- 2 * ext_encounter
    comment(new) <- "test"
    ext_encounter(params) <- new
    expect_equal(ext_encounter(params), new, ignore_attr = TRUE)
    expect_identical(comment(params@ext_encounter), "test")
})
