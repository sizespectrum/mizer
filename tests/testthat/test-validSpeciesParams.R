test_that("validSpeciesParams() works", {
    species_params <- NS_species_params
    
    # test w_mat
    sp <- species_params
    sp$w_mat[1] <- NA
    expect_message(sp <- validSpeciesParams(sp), NA)
    expect_equal(sp$w_mat[1], sp$w_max[1] / 4)
    sp$w_mat[2:4] <- 100
    expect_warning(sp <- validSpeciesParams(sp),
                   "For the species Sandeel, N.pout the value")
    expect_equal(sp$w_mat[2], sp$w_max[2] / 4)
    
    # test w_mat25
    sp <- species_params
    sp$w_mat25 <- c(NA, 1:11)
    expect_warning(validSpeciesParams(sp), NA)
    sp$w_mat25[2:5] <- 21
    expect_warning(sp <- validSpeciesParams(sp),
                   "For the species Sandeel, Dab the value")
    expect_true(is.na(sp$w_mat25[[2]]))
    expect_identical(sp$w_mat25[[3]], 21)
    
    # test w_min
    sp <- species_params
    sp$w_min <- c(NA, 1:11)
    expect_warning(validSpeciesParams(sp), NA)
    sp$w_min[2:5] <- 21
    expect_warning(sp <- validSpeciesParams(sp),
                   "For the species Sandeel, Dab the value")
    expect_identical(sp$w_min[[2]], 0.001)
    expect_identical(sp$w_min[[3]], 21)
    
    # test misspelling detection
    sp$wmat <- 1
    expect_warning(validSpeciesParams(sp),
                   "very close to standard parameter names") %>%
        suppressWarnings()
    
    # minimal species_params
    sp <- data.frame(species = c(2, 1),
                     w_max = c(100, 1000),
                     stringsAsFactors = FALSE)
    expect_s3_class(sp <- validSpeciesParams(sp), "data.frame")
    expect_equal(sp$w_mat, sp$w_max / 4)
    expect_equal(sp$alpha, c(0.6, 0.6))
    expect_equal(sp$interaction_resource, c(1, 1))
    expect_identical(rownames(sp), c("2", "1"))
})

test_that("validSpeciesParams converts from length to weight", {
    sp <- data.frame(species = 1:2,
                     l_max = 1:2,
                     a = 0.01, b = 3)
    sp2 <- validSpeciesParams(sp)
    expect_identical(sp2$w_max, c(0.01, 0.08))
})
