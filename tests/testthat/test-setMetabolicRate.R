params <- NS_params
no_sp <- nrow(params@species_params)

## setMetabolicRate ----
test_that("setMetabolicRate works", {
    expect_identical(setMetabolicRate(params, params@metab, 
                                      comment_metab = NULL), params)
    params@species_params$ks <- 2 * params@species_params$ks
    p2 <- setMetabolicRate(params)
    expect_identical(2 * params@metab, p2@metab)
})
test_that("setMetabolicRate can set exponent p", {
    # no change where p is already set in species_params
    params <- setMetabolicRate(params, p = 1)
    expect_identical(params@species_params$p, rep(0.7, 12))
    # but change where it is not
    params@species_params$p[[1]] <- NA
    params <- setMetabolicRate(params, p = 1)
    expect_identical(params@species_params$p, c(1, rep(0.7, 11)))
})
test_that("Comment works on metab", {
    metab <- params@metab
    comment(metab) <- "test"
    params <- setMetabolicRate(params, metab = metab)
    expect_identical(comment(params@metab), "test")
    
    # no message when nothing changes
    expect_message(setMetabolicRate(params), NA)
    # but message when a change is not stored due to comment
    params@species_params$k <- 1
    expect_message(setMetabolicRate(params),
                   "has been commented")
    
    # comment argument is ignored when there is a comment on intake_max
    params <- setMetabolicRate(params, metab = metab,
                               comment_metab = "overwrite")
    expect_identical(comment(params@metab), "test")
    # but it is used otherwise
    comment(metab) <- NULL
    params <- setMetabolicRate(params, metab = metab,
                               comment_metab = "overwrite")
    expect_identical(comment(params@metab), "overwrite")
})

# getMetabolicRate ----
test_that("getMetabolicRate works", {
    p <- setMetabolicRate(params, metab = getMetabolicRate(params), 
                          comment_metab = NULL)
    expect_identical(params, p)
})