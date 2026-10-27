# setBevertonHolt ----
NS_required_rdd <- getRequiredRDD(NS_params_small)
NS_rdd <- getRDD(NS_params_small)

test_that("setBevertonHolt sets erepro correctly when setting all values", {
    params <- setBevertonHolt(NS_params_small,
                              erepro = 10 * NS_params_small@species_params$erepro)
    expect_identical(params@species_params$erepro,
                     10 * NS_params_small@species_params$erepro)
    rdd <- getRDD(params)
    expect_equal(getRequiredRDD(params), rdd)
})

test_that("setBevertonHolt sets erepro correctly when setting same value for all species", {
    sp_name <- NS_params_small@species_params$species[3]
    expect_warning(params <- setBevertonHolt(NS_params_small, erepro = 0.1),
                   "For the following species `erepro` has been")
    expect_identical(params@species_params$R_max[params@species_params$species == sp_name],
                     Inf, ignore_attr = TRUE)
    expect_equal(NS_required_rdd, getRDD(params))
})

test_that("setBevertonHolt sets erepro correctly when setting same value for single species by using named vectors", {
    params <- NS_params_small
    erepro_old <- params@species_params$erepro
    erepro_new <- 0.1
    names(erepro_new) <- params@species_params$species[1]
    params <- setBevertonHolt(NS_params_small, erepro = erepro_new)
    expect_identical(params@species_params$erepro[2:3], erepro_old[2:3])
    expect_identical(params@species_params$erepro[1], 0.1, ignore_attr = TRUE)
})

test_that("setBevertonHolt sets erepro correctly when setting same value for some species by using named vector", {
    params <- NS_params_small
    erepro_old <- params@species_params$erepro
    erepro_new <- c(0.1, 0.2, NA)
    names(erepro_new) <- params@species_params$species[1:3]
    params <- setBevertonHolt(NS_params_small, erepro = erepro_new)
    expect_identical(params@species_params$erepro[3], erepro_old[3])
    expect_identical(params@species_params$erepro[1:2], c(0.1, 0.2), ignore_attr = TRUE)
    expect_equal(NS_required_rdd[1:2], getRDD(params)[1:2])
})

test_that("setBevertonHolt sets erepro correctly when setting same value for some species some species by NAs", {
    params <- NS_params_small
    erepro_old <- params@species_params$erepro
    erepro_new <- erepro_old
    erepro_new[3] <- NA
    erepro_new[1:2] <- c(0.1, 0.2)
    params <- setBevertonHolt(NS_params_small, erepro = erepro_new)
    expect_identical(params@species_params$erepro[3], erepro_old[3])
    expect_identical(params@species_params$erepro[1:2], c(0.1, 0.2), ignore_attr = TRUE)
})

test_that("setBevertonHolt updates `time_modified`", {
    expect_false(identical(setBevertonHolt(NS_params_small, erepro = 1)@time_modified,
                           NS_params_small@time_modified))
})

test_that("setBevertonHolt without extra arguments keeps previous erepro", {
    params <- setBevertonHolt(NS_params_small)
    expect_identical(params@species_params$erepro,
                     NS_params_small@species_params$erepro)
})

# R_max ----
test_that("setBevertonHolt sets R_max correctly when setting all values", {
    params <- setBevertonHolt(NS_params_small,
                              R_max = 10 * NS_params_small@species_params$R_max)
    expect_identical(params@species_params$R_max,
                     10 * NS_params_small@species_params$R_max)
    rdd <- getRDD(params)
    expect_equal(getRequiredRDD(params), rdd)
})

test_that("setBevertonHolt sets R_max only for named species when R_max is absent from given_species_params", {
    p <- NS_params
    given_species_params(p)$R_max <- NULL
    p2 <- setBevertonHolt(p, R_max = c(Sprat = 2e12))
    expect_equal(given_species_params(p2)$R_max[1], 2e12, ignore_attr = TRUE)
    expect_true(all(is.na(given_species_params(p2)$R_max[2:12])))
})

test_that("setBevertonHolt sets R_max correctly when setting values for some species by using named vector", {
    params <- NS_params_small
    R_max_old <- params@species_params$R_max
    R_max_new <- c(1e12, 1e13, NA)
    names(R_max_new) <- params@species_params$species[1:3]
    params <- setBevertonHolt(NS_params_small, R_max = R_max_new)
    expect_identical(params@species_params$R_max[3], R_max_old[3])
    expect_identical(params@species_params$R_max[1:2], c(1e12, 1e13), ignore_attr = TRUE)
})

test_that("setBevertonHolt issues warning when an R_max leads to an erepro > 1", {
    sp_name <- NS_params_small@species_params$species[3]
    rdi_sp <- getRDI(NS_params_small)[sp_name]
    req_rdd_sp <- getRequiredRDD(NS_params_small)[sp_name]
    # Set R_max just below the threshold where erepro > 1 is needed
    # Use a named single-species vector to avoid floating-point issues with other species
    threshold <- req_rdd_sp * rdi_sp / (rdi_sp - req_rdd_sp)
    R_max_new <- setNames(threshold * 0.99, sp_name)
    expect_warning(params <- setBevertonHolt(NS_params_small, R_max = R_max_new),
                   paste0("The following species require an unrealistic value greater than 1 for `erepro`: ", sp_name))
    expect_gt(params@species_params$erepro[params@species_params$species == sp_name], 1)
    expect_equal(params@species_params$R_max[3], unname(R_max_new), ignore_attr = TRUE)
})

# reproduction_level ----
test_that("setBevertonHolt and reproduction_level<- set reproduction_level correctly", {
    sp_name <- NS_params_small@species_params$species[3]
    expect_warning(params <- setBevertonHolt(NS_params_small, reproduction_level = 0.8),
                   paste0("The following species require an unrealistic value greater than 1 for `erepro`: ", sp_name))
    rdd <- getRDD(params)
    expect_equal(rdd, params@species_params$R_max * 0.8, ignore_attr = TRUE)
    expect_equal(getRequiredRDD(params), rdd)
    expect_equal(reproduction_level(params)[[1]], 0.8)

    # Test reproduction_level<-
    params2 <- NS_params_small
    expect_warning(reproduction_level(params2) <- 0.8,
                   paste0("The following species require an unrealistic value greater than 1 for `erepro`: ", sp_name))
    expect_equal(species_params(params2)$erepro, species_params(params)$erepro)
    expect_equal(species_params(params2)$R_max, species_params(params)$R_max)
})

test_that("reproduction_level is getRDD divided by R_max", {
    expect_equal(reproduction_level(NS_params_small),
                 NS_rdd / NS_params_small@species_params$R_max)
})

# R_factor ----
test_that("setBevertonHolt sets R_factor correctly", {
    sp_name <- NS_params_small@species_params$species[3]
    expect_warning(params <- setBevertonHolt(NS_params_small, R_factor = 1.5),
                   paste0("The following species require an unrealistic value greater than 1 for `erepro`: ", sp_name))
    rdd <- getRDD(params)
    expect_equal(rdd, params@species_params$R_max / 1.5, ignore_attr = TRUE)
    expect_equal(getRequiredRDD(params), rdd)
})

test_that("setRmax is an exact alias for setBevertonHolt", {
    p1 <- suppressWarnings(setRmax(NS_params_small, reproduction_level = 0.3))
    p2 <- suppressWarnings(setBevertonHolt(NS_params_small, reproduction_level = 0.3))
    expect_equal(species_params(p1)$erepro, species_params(p2)$erepro)
    expect_equal(species_params(p1)$R_max, species_params(p2)$R_max)
    expect_identical(reproduction_level(p1), reproduction_level(p2))
})

# special values ----
test_that("setBevertonHolt does nothing when called with only NA values", {
    erepro_new <- rep(NA, 3)
    params <- setBevertonHolt(NS_params_small, erepro = erepro_new)
    expect_identical(params, NS_params_small)
    erepro_new <- NA
    sp_name <- NS_params_small@species_params$species[3]
    names(erepro_new) <- sp_name
    params <- setBevertonHolt(NS_params_small, erepro = erepro_new)
    expect_identical(params, NS_params_small)
})

test_that("setBevertonHolt throws error when called with invalid values", {
    expect_error(setBevertonHolt(NS_params_small, erepro = rep("x", 3)),
                 "You provided invalid non-numeric values.")
    expect_error(setBevertonHolt(NS_params_small, erepro = 1, R_max = 2),
                 "You should only provide `params` and one other argument.")
    expect_error(setBevertonHolt(NS_params_small, R_max = 1:4),
                 "You need to supply a vector of length 3 or a single number or a named vector.")
    expect_error(setBevertonHolt(NS_params_small, reproduction_level = 2),
                 "The reproduction level must be smaller than 1 and non-negative.")
    expect_error(setBevertonHolt(NS_params_small, R_factor = 1/2),
                 "The R_factor must be greater than 1.")
})

test_that("reproduction_level of 0 works", {
    params <- setBevertonHolt(NS_params_small, reproduction_level = 0)
    expect_equal(getRDD(params), getRequiredRDD(params))
    expect_equal(params@species_params$R_max[1], Inf, ignore_attr = TRUE)

    params <- setBevertonHolt(NS_params_small, R_factor = Inf)
    expect_equal(getRDD(params), getRequiredRDD(params))
    expect_equal(params@species_params$R_max[1], Inf, ignore_attr = TRUE)
})

test_that("R_max is increased when needed", {
    sp_name1 <- NS_params_small@species_params$species[1]
    sp_name2 <- NS_params_small@species_params$species[2]
    expect_warning(p <- setBevertonHolt(NS_params_small, R_max = c(1, 2, NA)),
                 paste0("has been increased to give a reproduction level of 0.99: ", sp_name1, ", ", sp_name2))
    # R_max was increased from the requested values (1 and 2) to rdd/0.99
    expect_gt(p@species_params$R_max[1], 1)
    expect_gt(p@species_params$R_max[2], 2)
})

# absent species ----
test_that("setBevertonHolt leaves the parameters of an absent species alone", {
    params <- NS_params_small
    initialN(params)["Sprat", ] <- 0
    p <- setBevertonHolt(params)
    expect_identical(p@species_params$erepro, params@species_params$erepro)
    expect_identical(p@species_params$R_max[1], params@species_params$R_max[1])
    expect_equal(getRDD(p), getRequiredRDD(p))

    expect_no_warning(p <- setBevertonHolt(params, reproduction_level = 0))
    expect_identical(p@species_params$R_max[1], params@species_params$R_max[1])
    expect_identical(p@species_params$R_max[2:3], c(Inf, Inf),
                     ignore_attr = TRUE)
})

test_that("setBevertonHolt sets a given erepro or R_max on an absent species", {
    params <- NS_params_small
    initialN(params)["Sprat", ] <- 0
    p <- setBevertonHolt(params, erepro = c(Sprat = 0.3, Herring = 0.5))
    expect_identical(p@species_params$erepro[1:2], c(0.3, 0.5),
                     ignore_attr = TRUE)
    expect_identical(p@given_species_params$erepro[1], 0.3, ignore_attr = TRUE)
    expect_identical(p@species_params$R_max[1], params@species_params$R_max[1])
    expect_equal(getRDD(p)[1:2], getRequiredRDD(p)[1:2])

    p <- setBevertonHolt(params, R_max = c(Sprat = 1e10))
    expect_identical(p@species_params$R_max[1], 1e10, ignore_attr = TRUE)
    expect_identical(p@given_species_params$R_max[1], 1e10, ignore_attr = TRUE)
    expect_identical(p@species_params$erepro, params@species_params$erepro)
})

test_that("setBevertonHolt cannot give an absent species a reproduction level", {
    params <- NS_params_small
    initialN(params)["Sprat", ] <- 0
    expect_warning(p <- setBevertonHolt(params,
                                        reproduction_level = c(Sprat = 0.5)),
                   "absent.*: Sprat")
    expect_identical(p@species_params$erepro, params@species_params$erepro)
    expect_identical(p@species_params$R_max, params@species_params$R_max)
    expect_warning(setBevertonHolt(params, R_factor = c(Sprat = 2)),
                   "absent.*: Sprat")
    # Invalid values are rejected for an absent species too
    expect_error(setBevertonHolt(params, reproduction_level = c(Sprat = 2)),
                 "The reproduction level must be smaller than 1")
})

test_that("setBevertonHolt names a present species that has no reproduction", {
    params <- NS_params_small
    params@species_params$erepro[1] <- 0
    expect_error(setBevertonHolt(params), "no reproduction.*: Sprat\\.")
})

# info_level ----
test_that("info_level controls warnings in setBevertonHolt", {
    # info_level = 0 silences warnings
    expect_silent(setBevertonHolt(NS_params_small, reproduction_level = 0.8,
                                 info_level = 0))
    expect_silent(setBevertonHolt(NS_params_small, erepro = 0.1,
                                 info_level = 0))
    expect_silent(setBevertonHolt(NS_params_small, R_max = c(1, 2, NA),
                                 info_level = 0))

    # global option mizer_info_level = 0 silences warnings
    withr::with_options(list(mizer_info_level = 0), {
        expect_silent(setBevertonHolt(NS_params_small, reproduction_level = 0.8))
        p <- NS_params_small
        expect_silent(reproduction_level(p) <- 0.8)
    })

    # with_info_level deduplicates warnings
    sp_name <- NS_params_small@species_params$species[3]
    expect_warning(
        with_info_level({
            p <- setBevertonHolt(NS_params_small, reproduction_level = 0.8)
            setBevertonHolt(p, reproduction_level = 0.8)
        }),
        paste0("The following species require an unrealistic value greater than 1 for `erepro`: ", sp_name)
    )
})



