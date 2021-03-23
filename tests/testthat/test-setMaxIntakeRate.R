params <- NS_params
no_sp <- nrow(params@species_params)

## setMaxIntakeRate ----
test_that("ssetMaxIntakeRate works", {
    expect_identical(setMaxIntakeRate(params, params@intake_max, 
                                      comment_intake_max = NULL), params)
    params@species_params$h <- 2 * params@species_params$h
    p2 <- setMaxIntakeRate(params)
    expect_identical(2 * params@intake_max, p2@intake_max)
})
test_that("Comment works on intake_max", {
    intake_max <- params@intake_max
    comment(intake_max) <- "test"
    params <- setMaxIntakeRate(params, intake_max = intake_max)
    expect_identical(comment(params@intake_max), "test")
    
    # no message when nothing changes
    expect_message(setMaxIntakeRate(params), NA)
    # but message when a change is not stored due to comment
    params@species_params$h <- 1
    expect_message(setMaxIntakeRate(params),
                   "has been commented")
    
    # comment argument is ignored when there is a comment on intake_max
    params <- setMaxIntakeRate(params, intake_max = intake_max,
                               comment_intake_max = "overwrite")
    expect_identical(comment(params@intake_max), "test")
    # but it is used otherwise
    comment(intake_max) <- NULL
    params <- setMaxIntakeRate(params, intake_max = intake_max,
                               comment_intake_max = "overwrite")
    expect_identical(comment(params@intake_max), "overwrite")
})

# getMaxIntakeRate ----
test_that("getMaxIntakeRate works", {
    p <- setMaxIntakeRate(params, intake_max = getMaxIntakeRate(params), 
                          comment_intake_max = NULL)
    expect_identical(params, p)
})