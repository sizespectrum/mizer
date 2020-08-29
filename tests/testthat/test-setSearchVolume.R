params <- NS_params
no_sp <- nrow(params@species_params)

## setSearchVolume ----
test_that("setSearchVolume works", {
    expect_equal(setSearchVolume(params, params@search_vol), params)
    params@species_params$gamma <- 2 * params@species_params$gamma
    p2 <- setSearchVolume(params)
    expect_equal(2 * params@search_vol, p2@search_vol)
})
test_that("Comment works on search volume", {
    comment(params@search_vol) <- "test"
    params <- setSearchVolume(params, search_vol = params@search_vol)
    expect_identical(comment(params@search_vol), "test")
    expect_message(setSearchVolume(params), NA)
    params@species_params$gamma <- 1
    expect_message(setSearchVolume(params), "has been commented")
})

# getSearchVolume ----
test_that("getSearchVolume works", {
    p <- setSearchVolume(params, search_vol = getSearchVolume(params))
    expect_identical(params, p)
})