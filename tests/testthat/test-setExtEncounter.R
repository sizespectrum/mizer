test_that("setExtEncounter works", {
    params <- NS_params
    
    # Without ext_encounter argument, should be set to 0
    p2 <- setExtEncounter(params)
    zero_encounter <- p2@mu_b # to get right dimnames
    zero_encounter[] <- 0
    expect_identical(zero_encounter, p2@ext_encounter)
    
    # supplying ext_encounter
    p2 <- setExtEncounter(params, 3 * params@mu_b)
    expect_equal(p2@ext_encounter, 3 * params@mu_b)
    
    # only mu_b changed
    p2@ext_encounter <- params@ext_encounter
    expect_unchanged(p2, params)
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
    expect_identical(getExtEncounter(NS_params),
                     NS_params@ext_encounter)
})

test_that("Can get and set slot", {
    params <- NS_params
    ext_encounter <- getExtEncounter(params)
    expect_identical(ext_encounter(params), ext_encounter)
    new <- 2 * ext_encounter
    comment(new) <- "test"
    ext_encounter(params) <- new
    expect_identical(ext_encounter(params), new)
})
