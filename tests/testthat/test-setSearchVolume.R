params <- NS_params
no_sp <- nrow(params@species_params)

## setSearchVolume ----
test_that("setSearchVolume works", {
    expect_equal(setSearchVolume(params, params@search_vol, 
                                 comment = NULL), params)
    params@species_params$gamma <- 2 * params@species_params$gamma
    p2 <- setSearchVolume(params)
    expect_equal(2 * params@search_vol, p2@search_vol)
})
test_that("Comment works on search volume", {
    search_vol <- params@search_vol
    comment(search_vol) <- "test"
    params <- setSearchVolume(params, search_vol = search_vol)
    expect_identical(comment(params@search_vol), "test")
    
    # no message when nothing changes
    expect_message(setSearchVolume(params), NA)
    # but message when a change is not stored due to comment
    params@species_params$gamma <- 1
    expect_message(setSearchVolume(params), "has been commented")
    
    # comment argument is ignored when there is a comment on search_vol
    params <- setSearchVolume(params, search_vol = search_vol,
                               comment = "overwrite")
    expect_identical(comment(params@search_vol), "test")
    # but it is used otherwise
    comment(search_vol) <- NULL
    params <- setSearchVolume(params, search_vol = search_vol,
                              comment = "overwrite")
    expect_identical(comment(params@search_vol), "overwrite")
})

# getSearchVolume ----
test_that("getSearchVolume works", {
    p <- setSearchVolume(params, search_vol = getSearchVolume(params), 
                         comment = NULL)
    expect_identical(params, p)
})