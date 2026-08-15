# Opt-in gate for tests of experimental features.
#
# The direct (Newton) steady-state solver and the limit-cycle finder are both
# marked `lifecycle::badge("experimental")`, and their test files are by far the
# slowest in the suite (~74 s and ~62 s, a third of the total). They are
# therefore skipped unless MIZER_TEST_EXPERIMENTAL is set, so that the everyday
# `devtools::test()` loop stays quick. Run them with, from the shell,
#
#     MIZER_TEST_EXPERIMENTAL=true Rscript -e 'devtools::test()'
#
# or, within an R session,
#
#     Sys.setenv(MIZER_TEST_EXPERIMENTAL = "true")
#     devtools::test(filter = "steadyNewton|getLimitCycleSim")
#
# The full R CMD check workflows set the variable, so a deliberate check still
# covers this code. Remove the gate once these functions are no longer
# experimental.
testing_experimental <- function() {
    isTRUE(as.logical(Sys.getenv("MIZER_TEST_EXPERIMENTAL", "false")))
}

skip_unless_experimental <- function() {
    skip_if(!testing_experimental(),
            "Set MIZER_TEST_EXPERIMENTAL=true to test experimental features.")
}

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
# Additional cached objects — shared across test files to avoid rebuilding.
# R's copy-on-modify semantics ensure tests that mutate a local copy do not
# affect the cached originals.
#
# These are always built from scratch rather than loaded from a saved copy, so
# that every run exercises the constructors and no fixture can go stale against
# a change to the MizerParams class or to a rate setter's defaults (the problem
# that upgradeParams() exists to solve for saved objects).
#
# Each is built lazily, on first use, because most are needed by only a handful
# of test files while the whole helper is sourced by every one of them: eagerly
# they cost ~1.5 s, which is paid once per worker process when the suite runs
# in parallel. `delayedAssign()` keeps that transparent — test files still refer
# to the bare name — so a worker only pays for the fixtures its files touch.
delayedAssign("NS_sim_small",
    suppressMessages(project(NS_params_small, t_max = 3, t_save = 1,
                             progress_bar = FALSE)))
delayedAssign("single_sp_params", suppressMessages(newSingleSpeciesParams()))
delayedAssign("trait_params_small", suppressMessages(newTraitParams()))
delayedAssign("trait_params_2sp", suppressMessages(newTraitParams(no_sp = 2)))
delayedAssign("community_params_small", suppressMessages(newCommunityParams()))
# 3-species model with default no_w (differs from NS_params_small which uses no_w=20)
delayedAssign("NS_params_default_small", suppressMessages(
    newMultispeciesParams(NS_species_params_gears_small, inter_small,
                          info_level = 0)
))
# Single-species (Cod) model
delayedAssign("NS_params_cod_small", suppressMessages(
    newMultispeciesParams(NS_species_params_gears_small[3, ], info_level = 0)
))
# NS_params_small settled onto its steady state, for the many tests that need
# *a* model at steady state rather than to test steady() itself. The tolerance
# is the tightest any consumer needs (test-rate_functions.R checks the flux
# gradient against the mortality loss to 1e-8), and converging once to 1e-10 is
# cheaper than the several looser calls it replaces. Tests of steady() itself
# belong in test-steady.R and must of course keep calling it.
delayedAssign("NS_params_steady_small", suppressMessages(
    steady(NS_params_small, tol = 1e-10, t_max = 500,
           progress_bar = FALSE, info_level = 0)
))

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
