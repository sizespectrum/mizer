# Create an example MizerParams object
example_params <- function() {
    sp <- NS_species_params[9:11, ]
    # Make egg sizes different
    sp$w_min <- c(1e-3, 1e-2, 1e-1)

    # length-weight parameters
    sp$a <- c(0.01, 0.02, 0.03)
    sp$b <- c(3, 3, 3)

    gp <- data.frame(
        gear = c("Otter trawl", "Bottom trawl", "Bottom trawl"),
        species = c(sp$species[3], sp$species[3], sp$species[1]),
        catchability = c(0.1, 0.2, 0.3),
        sel_func = c("sigmoid_length", "knife_edge", "double_sigmoid_length"),
        knife_edge_size = c(NA, 40, NA),
        l50 = c(15, NA, 20),
        l25 = c(10, NA, 16),
        l50_right = c(NA, NA, 25),
        l25_right = c(NA, NA, 30)
    )

    params <- newMultispeciesParams(sp, gear_params = gp) |>
        suppressMessages()

    # Give diffusion to one species
    n <- params@species_params$n[1]
    d <- 0.1 * params@w^(n + 1)
    ext_diffusion(params)[1, ] <- d
    params
}

# Test that a MizerParams or MizerSim object has not changed except for the
# time_modified and perhaps a reordering of the species_params columns.
expect_unchanged <- function(object, expected) {
    if (is(object, "MizerParams")) {
        # has updated time_modified
        expect_false(identical(object@time_modified, expected@time_modified))
        object@time_modified <- expected@time_modified
        sp <- object@species_params
        sp_expected <- expected@species_params
    }
    if (is(object, "MizerSim")) {
        # has updated time_modified
        expect_false(identical(object@params@time_modified,
                               expected@params@time_modified))
        object@params@time_modified <- expected@params@time_modified
        sp <- object@params@species_params
        sp_expected <- expected@params@species_params
    }
    # Check that the species_params are unchanged except for a
    # reordering of the dataframe columns
    expect_equal(sp[, sort(names(sp))],
                     sp_expected[, sort(names(sp_expected))])
    # And if they were the same, remove the reordering
    if (is(object, "MizerParams")) {
        object@species_params <- sp_expected
    } else if (is(object, "MizerSim")) {
        object@params@species_params <- sp_expected
    }

    expect_equal(object, expected)
}
