local_edition(3)

# setMaxIntakeRate ----
test_that("setMaxIntakeRate works", {
    params <- NS_params
    params@species_params$h <- 2 * params@species_params$h
    p2 <- setMaxIntakeRate(params)
    expect_identical(2 * params@intake_max, p2@intake_max)
    # only intake max changed
    p2@intake_max <- params@intake_max
    p2@time_modified <- params@time_modified
    expect_identical(p2, params)
})
test_that("Comment works on intake_max", {
    params <- NS_params
    # if no comment, it is set automatically
    intake_max <- params@intake_max
    params <- setMaxIntakeRate(params, intake_max = intake_max)
    expect_identical(comment(params@intake_max), "set manually")
    
    # comment is stored
    comment(intake_max) <- "test"
    params <- setMaxIntakeRate(params, intake_max = intake_max)
    expect_identical(comment(params@intake_max), "test")
    
    # if no comment, previous comment is kept
    comment(intake_max) <- NULL
    params <- setMaxIntakeRate(params, intake_max = intake_max)
    expect_identical(comment(params@intake_max), "test")
    
    # no message when nothing changes
    expect_message(setMaxIntakeRate(params), NA)
    # but message when a change is not stored due to comment
    params@species_params$h <- 1
    expect_message(setMaxIntakeRate(params),  "has been commented")
    # Can reset
    p <- setMaxIntakeRate(params, reset = TRUE)
    expect_equal(p@intake_max[, 1], params@w[1]^params@species_params[["n"]],
                 ignore_attr = TRUE)
    expect_warning(setMaxIntakeRate(params, intake_max = intake_max,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

# getMaxIntakeRate ----
test_that("getMaxIntakeRate works", {
    expect_identical(getMaxIntakeRate(NS_params),
                     NS_params@intake_max)
})

test_that("Can get and set slot", {
    params <- NS_params
    intake_max <- getMaxIntakeRate(params)
    expect_identical(intake_max(params), intake_max)
    new <- 2 * intake_max
    comment(new) <- "test"
    intake_max(params) <- new
    expect_identical(intake_max(params), new)
})
