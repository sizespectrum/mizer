# setSearchVolume ----
test_that("setSearchVolume works", {
    params <- NS_params
    params@species_params$gamma <- 2 * params@species_params$gamma
    p2 <- setSearchVolume(params)
    expect_equal(2 * params@search_vol, p2@search_vol)
    # only search_vol changed
    p2@search_vol <- params@search_vol
    expect_unchanged(p2, params)
})
test_that("Comment works on search volume", {
    params <- NS_params
    # if no comment, it is set automatically
    search_vol <- params@search_vol
    params <- setSearchVolume(params, search_vol = search_vol)
    expect_identical(comment(params@search_vol), "set manually")
    
    # comment is stored
    comment(search_vol) <- "test"
    params <- setSearchVolume(params, search_vol = search_vol)
    expect_identical(comment(params@search_vol), "test")
    
    # if no comment, previous comment is kept
    comment(search_vol) <- NULL
    params <- setSearchVolume(params, search_vol = search_vol)
    expect_identical(comment(params@search_vol), "test")
    
    # no message when nothing changes
    expect_message(setSearchVolume(params), NA)
    # but message when a change is not stored due to comment
    params@species_params$gamma <- 1
    expect_message(setSearchVolume(params),  "has been commented")
    # Can reset
    p <- setSearchVolume(params, reset = TRUE)
    expect_equal(p@search_vol[, 1], params@w[1]^params@species_params$q,
                 check.attributes = FALSE)
    expect_warning(setSearchVolume(params, search_vol = search_vol,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

# getSearchVolume ----
test_that("getSearchVolume works", {
    expect_identical(getSearchVolume(NS_params),
                     NS_params@search_vol)
})


test_that("Can get and set search_vol slot", {
    params <- NS_params
    new <- 2 * search_vol(params)
    comment(new) <- "test"
    search_vol(params) <- new
    expect_identical(search_vol(params), new)
})
