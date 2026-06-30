# Create an example MizerParams object
example_params <- function() {
    sp <- NS_species_params_small
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

# Build small 3-species (Sprat, Herring, Cod), 3-gear, 20-bin versions for faster tests.
# These are named _small to avoid overwriting the package datasets in the global environment.
# Load package datasets into a temp env to avoid polluting the test environment.
.mizer_test_data <- new.env(parent = emptyenv())
data("NS_species_params", package = "mizer", envir = .mizer_test_data)
data("NS_species_params_gears", package = "mizer", envir = .mizer_test_data)
data("inter", package = "mizer", envir = .mizer_test_data)
NS_species_params_small <- .mizer_test_data$NS_species_params[c(1, 4, 11), ]
NS_species_params_gears_small <- .mizer_test_data$NS_species_params_gears[
    .mizer_test_data$NS_species_params_gears$species %in% NS_species_params_small$species, ]
inter_small <- .mizer_test_data$inter[
    NS_species_params_small$species, NS_species_params_small$species]
rm(.mizer_test_data)
NS_params_small <- suppressMessages(
    newMultispeciesParams(NS_species_params_gears_small, inter_small, no_w = 20, info_level = 0)
)
# Mirror the given_species_params that calibration sets in the package NS_params.
NS_params_small@given_species_params$gamma <- NS_params_small@species_params$gamma
NS_params_small@given_species_params$f0 <- rep(0.6, nrow(NS_params_small@species_params))
NS_params_small@species_params$f0 <- NS_params_small@given_species_params$f0
NS_params_small@given_species_params$h <- NS_params_small@species_params$h
NS_params_small@given_species_params$ks <- NS_params_small@species_params$ks
# Set non-zero initial effort (matches pattern in original NS_params)
initial_effort(NS_params_small) <- c(Industrial = 0, Pelagic = 1, Otter = 0.5)
# Create NS_sim_small to match the 3-species NS_params_small
NS_sim_small <- suppressMessages(project(NS_params_small, t_max = 3, t_save = 1, progress_bar = FALSE))

# Additional cached objects — shared across test files to avoid rebuilding.
# R's copy-on-modify semantics ensure tests that mutate a local copy do not
# affect the cached originals.
single_sp_params <- suppressMessages(newSingleSpeciesParams())
trait_params_small <- suppressMessages(newTraitParams())
trait_params_2sp <- suppressMessages(newTraitParams(no_sp = 2))
community_params_small <- suppressMessages(newCommunityParams())
# 3-species model with default no_w (differs from NS_params_small which uses no_w=20)
NS_params_default_small <- suppressMessages(
    newMultispeciesParams(NS_species_params_gears_small, inter_small, info_level = 0)
)
# Single-species (Cod) model
NS_params_cod_small <- suppressMessages(
    newMultispeciesParams(NS_species_params_gears_small[3, ], info_level = 0)
)

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
    # reordering of the dataframe columns (ignore column-level name attributes
    # since named numeric columns are new behaviour from $.species_params)
    expect_equal(sp[, sort(names(sp))],
                     sp_expected[, sort(names(sp_expected))],
                     ignore_attr = TRUE)
    # And if they were the same, remove the reordering; also align class of
    # given_species_params so old fixtures with plain data.frame compare equal
    if (is(object, "MizerParams")) {
        object@species_params <- sp_expected
        object@given_species_params <- expected@given_species_params
    } else if (is(object, "MizerSim")) {
        object@params@species_params <- sp_expected
        object@params@given_species_params <- expected@params@given_species_params
    }

    expect_equal(object, expected, ignore_attr = TRUE)
}
