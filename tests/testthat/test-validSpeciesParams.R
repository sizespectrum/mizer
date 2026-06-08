test_that("validSpeciesParams() works", {
    species_params <- NS_species_params
    
    # test w_mat
    sp <- species_params
    sp$w_mat[1] <- NA
    expect_message(sp <- validSpeciesParams(sp), NA)
    expect_equal(sp$w_mat[1], sp$w_max[1] / 4)
    # Set w_mat > w_max for species 1 and 2 (Sprat w_max=33, Herring w_max=334)
    sp$w_mat[1:2] <- sp$w_max[1:2] * 2
    sp1 <- sp$species[1]; sp2 <- sp$species[2]
    expect_warning(sp <- validSpeciesParams(sp),
                   paste0("For the species ", sp1, ", ", sp2, " the value"))
    expect_equal(sp$w_mat[1], sp$w_max[1] / 4)

    # test w_mat25
    sp <- species_params
    sp$w_mat25 <- c(NA, 1, 2)
    expect_warning(validSpeciesParams(sp), NA)
    # Set w_mat25 > w_mat for species 1 and 2
    sp$w_mat25[1:2] <- sp$w_mat[1:2] * 2
    expect_warning(sp <- validSpeciesParams(sp),
                   paste0("For the species ", sp1, ", ", sp2, " the value"))
    expect_true(is.na(sp$w_mat25[[1]]))
    expect_identical(sp$w_mat25[[3]], 2)

    # test w_min
    sp <- species_params
    sp$w_min <- c(NA, 1, 2)
    expect_warning(validSpeciesParams(sp), NA)
    # Set w_min > w_max for species 1 and 2
    sp$w_min[1:2] <- sp$w_max[1:2] * 2
    expect_warning(sp <- validSpeciesParams(sp),
                   paste0("For the species ", sp1, ", ", sp2, " the value"))
    expect_identical(sp$w_min[[1]], 0.001)
    expect_identical(sp$w_min[[3]], 2)
    
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

test_that("validGivenSpeciesParams checks documented error cases and signals inconsistency info", {
    expect_error(
        validGivenSpeciesParams(data.frame(w_max = 1)),
        "needs a column 'species'"
    )
    expect_error(
        validGivenSpeciesParams(data.frame(species = c("a", "a"), w_max = c(1, 2))),
        "multiple rows for the same species"
    )
    expect_error(
        validGivenSpeciesParams(data.frame(species = c("a", "b"), w_max = c(1, NA))),
        "specify maximum sizes for all species"
    )

    sp <- data.frame(species = c("a", "b"),
                     l_max = c(10, 20),
                     w_max = c(1, 1000),
                     a = 0.01,
                     b = 3)
    expect_condition(
        validGivenSpeciesParams(sp),
        "For the following species I will ignore your value for l_max",
        class = "info_about_default"
    )
})

test_that("validSpeciesParams sets the documented defaults", {
    sp <- data.frame(species = c("a", "b"),
                     w_max = c(10, 100),
                     w_mat = c(NA, 20),
                     w_min = c(NA, 1),
                     alpha = c(NA, 0.5),
                     interaction_resource = c(NA, 2),
                     n = c(NA, 0.8),
                     p = c(NA, 0.7))

    sp2 <- validSpeciesParams(sp)

    expect_equal(sp2$w_repro_max, sp2$w_max)
    expect_equal(sp2$w_mat, c(10 / 4, 20))
    expect_equal(sp2$w_min, c(0.001, 1))
    expect_equal(sp2$alpha, c(0.6, 0.5))
    expect_equal(sp2$interaction_resource, c(1, 2))
    expect_equal(sp2$n, c(3/4, 0.8))
    expect_equal(sp2$p, c(3/4, 0.7))
})

test_that("completeSpeciesParams is a deprecated alias for validSpeciesParams", {
    sp <- data.frame(species = c("a", "b"),
                     w_max = c(10, 100),
                     stringsAsFactors = FALSE)

    expect_identical(completeSpeciesParams(sp), validSpeciesParams(sp))
})
