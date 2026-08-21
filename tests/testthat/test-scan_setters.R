# Setters --------------------------------------------------------------------

test_that("scanEffort() sets the effort of the right gears", {
    p <- scanEffort()(NS_params_small, 2)
    expect_true(all(initial_effort(p) == 2))

    p <- scanEffort("Otter")(NS_params_small, 2)
    expect_identical(unname(initial_effort(p)[["Otter"]]), 2)
    expect_identical(initial_effort(p)[c("Industrial", "Pelagic")],
                     initial_effort(NS_params_small)[c("Industrial", "Pelagic")])

    expect_error(scanEffort("Nonexistent")(NS_params_small, 1),
                 "no gear called Nonexistent")
})

test_that("scanFishingMortality() is idempotent", {
    # This is what lets it be used with continuation, where it is applied to the
    # object it returned at the previous scan value.
    f <- scanFishingMortality("Cod")
    p1 <- suppressMessages(f(NS_params_small, 0.3))
    p2 <- suppressMessages(f(p1, 0.3))

    expect_equal(gear_params(p2), gear_params(p1))
    expect_equal(initial_effort(p2), initial_effort(p1))
    # Only the fishing mortality values matter here; the params attribute that
    # getFMort() carries differs in its time_modified stamp.
    expect_equal(unclass(getFMort(p2))[, ], unclass(getFMort(p1))[, ],
                 ignore_attr = TRUE)
    # The extra gear is installed exactly once.
    added <- setdiff(dimnames(p2@catchability)[[1]],
                     dimnames(NS_params_small@catchability)[[1]])
    expect_length(added, 1)
})

test_that("a setter reused on another model does not hijack its gear", {
    # The setter carries the name of the gear it installed. Handed a different
    # model that happens to have a gear of that name catching the same species
    # with catchability 1, it must not take that for its own installation: the
    # original gears would then still be fishing and the scanned mortality
    # would add to them instead of replacing them.
    f <- scanFishingMortality("Cod")
    p_a <- suppressMessages(f(NS_params_small, 0.3))
    installed <- setdiff(mizer:::gear_names(p_a),
                         mizer:::gear_names(NS_params_small))
    expect_length(installed, 1)

    # A second model carrying a decoy gear of exactly that name
    p_b <- NS_params_small
    gp <- gear_params(p_b)
    decoy <- gp[gp$species == "Cod", ][1, ]
    decoy$gear <- installed
    decoy$catchability <- 1
    gear_params(p_b) <- rbind(gp, decoy)
    initial_effort(p_b)[installed] <- 1

    # Applied to that model the setter must still deliver the mortality asked
    # for: doubling the request doubles Cod's fishing mortality.
    r1 <- suppressMessages(f(p_b, 0.3))
    r2 <- suppressMessages(f(p_b, 0.6))
    expect_equal(max(getFMort(r2)["Cod", ]), 2 * max(getFMort(r1)["Cod", ]))
    # and other species are untouched
    expect_equal(unclass(getFMort(r2))["Herring", ],
                 unclass(getFMort(r1))["Herring", ])
})

test_that("a reused setter still requires its requested gear", {
    f <- scanFishingMortality("Cod", gear = "Otter")
    p_a <- suppressMessages(f(NS_params_small, 0.3))
    installed <- setdiff(mizer:::gear_names(p_a),
                         mizer:::gear_names(NS_params_small))
    expect_length(installed, 1)

    # This model has a decoy bearing the remembered scan-gear name, but the
    # original Cod--Otter row that proves a genuine installation is absent.
    p_b <- NS_params_small
    gp <- gear_params(p_b)
    gp <- gp[gp$species != "Cod" | gp$gear != "Otter", ]
    decoy <- gear_params(p_a)
    decoy <- decoy[decoy$gear == installed, ]
    gear_params(p_b) <- rbind(gp, decoy)
    initial_effort(p_b)[installed] <- 1

    expect_error(suppressMessages(f(p_b, 0.6)),
                 "The gear Otter does not catch Cod")
})

test_that("scanFishingMortality() does not hijack a gear of its own name", {
    # A model that already has a gear called "scan" must not have that gear's
    # effort scanned in place of the target species' fishing mortality.
    p <- NS_params_small
    gp <- gear_params(p)
    decoy <- gp[gp$species == "Herring", ][1, ]
    decoy$gear <- "scan"
    gear_params(p) <- rbind(gp, decoy)
    initial_effort(p)["scan"] <- 1

    f <- scanFishingMortality("Cod")
    p1 <- suppressMessages(f(p, 0.3))
    p2 <- suppressMessages(f(p1, 0.6))

    # Doubling the requested value doubles Cod's fishing mortality ...
    expect_equal(max(getFMort(p2)["Cod", ]), 2 * max(getFMort(p1)["Cod", ]))
    # ... and leaves the pre-existing "scan" gear's effort alone.
    expect_identical(initial_effort(p2)[["scan"]],
                     initial_effort(p)[["scan"]])
})

test_that("scanFishingMortality() varies only the target species", {
    f <- scanFishingMortality("Cod")
    p <- suppressMessages(f(NS_params_small, 0.3))
    low <- getFMort(f(p, 0.3))
    high <- getFMort(f(p, 0.6))
    expect_equal(max(high["Cod", ]), 2 * max(low["Cod", ]))
    expect_equal(unclass(high)["Herring", ], unclass(low)["Herring", ])
})

test_that("scanSpeciesParam() sets a parameter and propagates the change", {
    # A parameter the user supplied, so present in given_species_params().
    # mizer warns that w_mat25 no longer sits below w_mat, which is the
    # validation now running as it should.
    p <- suppressWarnings(suppressMessages(
        scanSpeciesParam("Cod", "w_mat")(NS_params_small, 100)))
    expect_identical(species_params(p)[["w_mat"]][[3]], 100)
    expect_identical(given_species_params(p)[["w_mat"]][[3]], 100)

    # A parameter mizer defaulted, so absent from given_species_params(). This
    # used to build a malformed column and then fail in the rate calculations.
    expect_false("alpha" %in% names(given_species_params(NS_params_small)))
    p <- suppressMessages(scanSpeciesParam("Cod", "alpha")(NS_params_small, 0.4))
    expect_identical(species_params(p)[["alpha"]][[3]], 0.4)
    expect_type(species_params(p)[["alpha"]], "double")
    # and it reaches the rates that depend on it
    expect_false(isTRUE(all.equal(unclass(getEGrowth(p))[3, ],
                                  unclass(getEGrowth(NS_params_small))[3, ])))

    # An unknown parameter is rejected rather than silently creating a column
    expect_error(scanSpeciesParam("Cod", "not_a_parameter")(NS_params_small, 1),
                 "no species parameter called")
    expect_error(scanSpeciesParam("Nonexistent", "w_mat")(NS_params_small, 1),
                 "no species called")
})
