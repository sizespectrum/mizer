
example_manipulate_params <- example_params()
small_trait_grid_manipulate_params <- newTraitParams(no_sp = 2, no_w = 36,
                                                     info_level = 0)
ns_manipulate_params <- newMultispeciesParams(NS_species_params_small,
                                              info_level = 0)
ns_manipulate_params_36 <- newMultispeciesParams(
    NS_species_params_small,
    no_w = 36,
    max_w = 39900,
    min_w_pp = 9e-14,
    info_level = 0
)
remove_ns_species <- NS_species_params_small$species[2]
reduced_ns_species_params <-
    NS_species_params_small[NS_species_params_small$species != remove_ns_species, ]
reduced_ns_manipulate_params_36 <- newMultispeciesParams(
    reduced_ns_species_params,
    no_w = 36,
    max_w = 39900,
    min_w_pp = 9e-14,
    info_level = 0
)
lower_ns_species_params <- NS_species_params_small
lower_ns_species_params$species <- tolower(lower_ns_species_params$species)
lower_ns_manipulate_params <- newMultispeciesParams(lower_ns_species_params,
                                                    info_level = 0)

# addSpecies ----
test_that("addSpecies works when adding a second identical species", {
    p <- trait_params_small
    no_sp <- nrow(p@species_params)
    species_params <- p@species_params[5, ]
    species_params$species <- "new"
    # Adding species 5 again should lead two copies of the species
    pa <- addSpecies(p, species_params)
    expect_identical(pa@metab[5, ], pa@metab[no_sp + 1, ])
    expect_identical(pa@psi[5, ], pa@psi[no_sp + 1, ])
    expect_identical(pa@ft_pred_kernel_e[5, ], pa@ft_pred_kernel_e[no_sp + 1, ])

    # test that we can remove species again
    pr <- removeSpecies(pa, "new")

})
test_that("addSpecies does not allow duplicate species", {
    p <- example_manipulate_params
    species_params <- p@species_params[3, ]
    expect_error(addSpecies(p, species_params),
                 "You can not add species that are already there.")
})
test_that("addSpecies handles gear params correctly", {
    p <- trait_params_2sp
    sp <- data.frame(species = c("new1", "new2"),
                     w_max = c(10, 100),
                     k_vb = c(4, 1),
                     n = 2 / 3,
                     p = 2 / 3)
    gp <- data.frame(gear = c("gear1", "gear2", "gear1"),
                     species = c("new1", "new2", "new2"),
                     sel_func = "knife_edge",
                     knife_edge_size = c(5, 5, 50))

    # If no initial_effort for new gear is provided, it is the edition default
    # Wrapping in `expect_warning()` to ignore warnings about unrealistic
    # reproductive efficiency
    (pa <- addSpecies(p, sp, gp)) |>
        expect_message() |>
        expect_warning()
    default_effort <- ifelse(defaults_edition() < 2, 0, 1)
    expect_identical(pa@initial_effort,
                     c(knife_edge_gear = default_effort,
                       gear1 = default_effort, gear2 = default_effort))
    expect_identical(nrow(pa@gear_params), 5L)

    # effort for existing gear is not changed
    extra_effort <- c(gear1 = 2, gear2 = 3)
    (pa <- addSpecies(p, sp, gp, initial_effort = extra_effort)) |>
        expect_message() |>
        expect_warning()
    expect_identical(pa@initial_effort, c(knife_edge_gear = default_effort, extra_effort))

    effort <- 2
    addSpecies(p, sp, gp, initial_effort = effort) |>
        expect_message() |>
        expect_error("The `initial_effort` must be a named list or vector")

    effort <- c(knife_edge_gear = 1)
    addSpecies(p, sp, gp, initial_effort = effort) |>
        expect_message() |>
        expect_error("The names of the `initial_effort` do not match the names of the new gears.")
})

test_that("addSpecies handles interaction matrix correctly", {
    p <- trait_params_2sp
    p <- setInteraction(p, interaction = matrix(1:4/8, ncol = 2))
    sp <- data.frame(species = c("new1", "new2"),
                     w_max = c(10, 100),
                     k_vb = c(4, 1),
                     n = 2/3,
                     p = 2/3)

    interaction <- matrix(1:4/4, ncol = 2)
    ones <- matrix(rep(1, 4), ncol = 2)
    (pa <- addSpecies(p, sp, interaction = interaction)) |>
        expect_message() |>
        expect_warning()
    expect_equal(pa@interaction[3:4, 3:4], interaction, ignore_attr = TRUE)
    expect_equal(pa@interaction[1:2, 3:4], ones, ignore_attr = TRUE)
    expect_equal(pa@interaction[3:4, 1:2], ones, ignore_attr = TRUE)
    expect_equal(pa@interaction[1:2, 1:2], p@interaction, ignore_attr = TRUE)

    interaction <- matrix(1:16/16, ncol = 4)
    (pa <- addSpecies(p, sp, interaction = interaction)) |>
        expect_message() |>
        expect_warning("The following species require an unrealistic value greater than 1 for `erepro`: new2")
    expect_equal(pa@interaction, interaction, ignore_attr = TRUE)

    addSpecies(p, sp, interaction = matrix(1:9, ncol = 3)) |>
        expect_warning() |>
        expect_error("Interaction matrix has invalid dimensions.")
})
test_that("addSpecies works when adding a species with a larger w_max", {
    sp <- data.frame(species = "Blue whale", w_max = 5e4,
                     w_mat = 1e3, beta = 1000, sigma = 2,
                     k_vb = 0.6, gear = 'Whale hunter')
    params <- example_manipulate_params
    # change a slot to test that such changes will be preserved
    params <- setMaxIntakeRate(params, 2 * getMaxIntakeRate(params))

    (p <- addSpecies(params, sp)) |>
        expect_message()
    expect_identical(p@w[1:100], params@w)
    expect_identical(p@w_full[seq_along(params@w_full)], params@w_full)
    expect_lte(5e4, max(p@w))
    # changed rates are preserved
    expect_equal(getMaxIntakeRate(p)[1:3, 1:100],
                 getMaxIntakeRate(params), ignore_attr = TRUE)
})
test_that("addSpecies works when adding a species with a smaller w_min", {
    sp <- data.frame(species = "Blue whale", w_max = 5e4, w_min = 1e-5,
                     w_mat = 1e3, beta = 1000, sigma = 2,
                     k_vb = 0.6, gear = 'Whale hunter')
    params <- NS_params_small
    # change a slot to test that such changes will be preserved
    params <- setMaxIntakeRate(params, 2 * getMaxIntakeRate(params))

    (p <- addSpecies(params, sp)) |>
        expect_message()
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    w_start <- which(abs(p@w - params@w[1]) < 1e-10)[1]
    w_end <- w_start + no_w - 1
    expect_equal(p@w[w_start:w_end], params@w)
    expect_equal(p@w_full[seq_along(params@w_full)], params@w_full)
    expect_gte(1e-5, min(p@w))
    # changed rates are preserved
    expect_equal(getMaxIntakeRate(p)[1:no_sp, w_start:w_end],
                 getMaxIntakeRate(params), ignore_attr = TRUE)
})

test_that("addSpecies has other documented properties", {
    sp <- data.frame(species = c("new1", "new2"),
                     w_max = c(10, 100),
                     k_vb = c(4, 1),
                     n = 2 / 3,
                     p = 2 / 3)
    (p <- addSpecies(example_manipulate_params, sp)) |>
        expect_message()

    # New species have 0 reproduction level
    expect_equal(getReproductionLevel(p)[4:5],
                 c(new1 = 1 / 4, new2 = 1 / 4))

    # Maximum of ratio between new species density and Sheldon density is 1/100
    fraction <- p@initial_n[4, ] /
        (p@resource_params$kappa * p@w ^ -p@resource_params$lambda)
    expect_equal(max(fraction), 1 / 100)
    expect_false(identical(p@time_modified, example_manipulate_params@time_modified))
})

test_that("Added species stay at low abundance", {
    # Use example from man page
    params <- trait_params_small
    species_params <- data.frame(
        species = "mullet",
        w_max = 173,
        w_mat = 15,
        beta = 283,
        sigma = 1.8,
        k_vb = 0.6,
        a = 0.0085,
        b = 3.11
    )
    expect_message(params <- addSpecies(params, species_params))
    no_sp <- nrow(params@species_params)
    sim <- project(params, t_max = 1, progress_bar = FALSE)
    expect_lt(finalN(sim)[no_sp, 1] / initialN(sim)[no_sp, 1], 1.05)
})

test_that("addSpecies preserves both given and other species params", {

    params <- trait_params_small
    params@given_species_params$b <- 3
    params@species_params$w_mat25 <- params@species_params$w_mat25 * 1.01
    sp <- data.frame(species = "new",
                     w_max = 10, test = "test")
    expect_message(p <- addSpecies(params, sp))
    no_sp <- nrow(params@species_params)
    expect_identical(p@species_params$w_mat25[1:no_sp],
                     params@species_params$w_mat25[1:no_sp])
    expect_identical(p@given_species_params$b[1:no_sp],
                     params@given_species_params$b[1:no_sp])
    expect_in("test", names(p@species_params))
    expect_in("test", names(p@given_species_params))
})

randomise_upgrade_preserved_slots <- function(params) {
    component <- "random_component"
    zero_rate <- function(params, n, n_pp, n_other, component, ...) {
        params@ext_encounter * 0
    }
    zero_mort <- function(params, n, n_pp, n_other, component, ...) {
        params@mu_b * 0
    }

    params@ext_diffusion[] <- runif(length(params@ext_diffusion), max = 1e-6)
    params@initial_n_other <- stats::setNames(list(runif(4)), component)
    params@other_dynamics <- stats::setNames(list("random_dynamics"), component)
    params@other_encounter <- stats::setNames(list(zero_rate), component)
    params@other_mort <- stats::setNames(list(zero_mort), component)
    params@other_params <- stats::setNames(list(list(values = runif(3))),
                                           component)
    params@other_params$other <- list(values = runif(2))
    params@use_predation_diffusion <- TRUE
    attr(params@rates_funcs, "random") <- runif(2)
    params
}

test_that("addSpecies preserves slots recently added to upgradeParams", {
    set.seed(42)
    params <- small_trait_grid_manipulate_params
    params <- randomise_upgrade_preserved_slots(params)

    sp <- data.frame(species = "new", w_max = 10, k_vb = 1)
    p <- suppressWarnings(suppressMessages(
        addSpecies(params, sp, info_level = 0)
    ))
    old_sp <- seq_len(nrow(params@species_params))
    old_w <- seq_along(params@w)

    expect_equal(p@ext_diffusion[old_sp, old_w], params@ext_diffusion,
                 ignore_attr = TRUE)
    for (slot in c("initial_n_other", "other_dynamics", "other_encounter",
                   "other_mort", "other_params", "rates_funcs",
                   "use_predation_diffusion")) {
        expect_identical(slot(p, slot), slot(params, slot))
    }
})

test_that("addSpecies preserves the class", {
    if (!methods::isClass("AddSpeciesTestParams")) {
        methods::setClass("AddSpeciesTestParams", contains = "MizerParams")
    }
    params <- as(small_trait_grid_manipulate_params, "AddSpeciesTestParams")
    sp <- data.frame(species = "new", w_max = 10, k_vb = 1)
    p <- suppressWarnings(suppressMessages(
        addSpecies(params, sp, info_level = 0)
    ))

    expect_s4_class(p, "AddSpeciesTestParams")
    expect_true(validObject(p))
})

# removeSpecies ----
test_that("removeSpecies works", {
    remove <- remove_ns_species
    params <- ns_manipulate_params_36
    p1 <- removeSpecies(params, species = remove)
    expect_equal(nrow(p1@species_params), nrow(params@species_params) - length(remove))
    p2 <- reduced_ns_manipulate_params_36
    p2@linecolour[2] = "#8da600" # update line colour
    expect_unchanged(p1, p2)
    sim1 <- project(p1, t_max = 0.4, t_save = 0.4)
    sim2 <- project(p2, t_max = 0.4, t_save = 0.4)
    expect_identical(sim1@n[2, 2, ], sim2@n[2, 2, ])
})
test_that("removeSpecies works with 3d pred kernel", {
    # It should make no difference whether we first set full pred kernel and
    # then remove a species, or the other way around.
    params1 <- example_manipulate_params
    sp_name <- params1@species_params$species[3]
    params1 <- setPredKernel(params1, pred_kernel = getPredKernel(params1))
    params1 <- removeSpecies(params1, sp_name)
    params2 <- example_manipulate_params
    params2 <- removeSpecies(params2, sp_name)
    params2 <- setPredKernel(params2, pred_kernel = getPredKernel(params2))
    expect_unchanged(params1, params2)
})
test_that("removeSpecies works correctly on gear_params", {
    p <- example_manipulate_params
    sp_name <- p@species_params$species[3]
    params <- removeSpecies(p, sp_name)
    expect_equal(nrow(params@gear_params), 1)
})

test_that("removeSpecies accepts numeric and logical selectors", {
    by_name <- removeSpecies(NS_params_small, c("Cod", "Herring"))
    by_index <- removeSpecies(NS_params_small, c(2, 3))
    by_logical <- removeSpecies(NS_params_small,
                                species_params(NS_params_small)$species %in%
                                    c("Cod", "Herring"))
    expect_equal(by_index, by_name, ignore_attr = TRUE)
    expect_equal(by_logical, by_name, ignore_attr = TRUE)
    expect_false(identical(by_name@time_modified, NS_params_small@time_modified))
})

test_that("adding and then removing species leaves params unaltered", {
    params <- ns_manipulate_params
    # two arbitrary species
    sp <- data.frame(species = c("new1", "new2"),
                     w_max = c(10, 100),
                     k_vb = c(4, 2),
                     test = "test",
                     stringsAsFactors = FALSE)
    # add comments to test that they will be preserved as well
    comment(params) <- "test"
    for (slot in (slotNames(params))) {
        comment(slot(params, slot)) <- slot
    }
    # But no comments in fields that would disable addSpecies
    comment(params@pred_kernel) <- NULL
    comment(params@catchability) <- NULL
    comment(params@selectivity) <- NULL
    (params2 <- addSpecies(params, sp) |>
        removeSpecies(c("new1", "new2"))) |>
        expect_message()

    # For now the linecolour and linetype are not preserved
    # TODO: fix this in the next overhaul of linecolour and linetype code
    params2@linecolour <- params@linecolour
    params2@linetype <- params@linetype
    params2@species_params$linecolour <- NULL
    params2@species_params$linetype <- NULL
    params2@given_species_params$linecolour <- NULL
    params2@given_species_params$linetype <- NULL
    # erepro may be added to given_species_params by setBevertonHolt during addSpecies
    params2@given_species_params$erepro <- NULL
    # comment on w_min_idx are not preserved
    comment(params@w_min_idx) <- NULL
    expect_unchanged(params, params2)
})

# renameSpecies ----
test_that("renameSpecies works", {
    sp <- lower_ns_species_params
    p <- ns_manipulate_params
    replace <- NS_species_params_small$species
    names(replace) <- sp$species
    p2 <- lower_ns_manipulate_params
    p2 <- renameSpecies(p2, replace)
    p2@time_modified <- p@time_modified
    p2@time_created <- p@time_created
    expect_identical(p, p2)
})

test_that("renameSpecies updates linked species names", {
    replace <- c(Cod = "Kabeljau", Herring = "Hering")
    p <- renameSpecies(NS_params_small, replace)
    expect_false(identical(p@time_modified, NS_params_small@time_modified))
    expect_true(all(replace %in% species_params(p)$species))
    expect_true(all(replace %in% names(getColours(p))))
    expect_true(all(replace %in% names(getLinetypes(p))))
    expect_true(all(replace %in% gear_params(p)$species))
    expect_true(all(replace %in% dimnames(initialN(p))$sp))
    expect_true(all(replace %in% dimnames(getSelectivity(p))$sp))
    expect_true(all(replace %in% dimnames(getCatchability(p))$sp))
})
test_that("renameSpecies warns on wrong names", {
    expect_error(renameSpecies(example_manipulate_params,
                               c(Kod = "cod", Hadok = "haddock")),
                 "Kod, Hadok do not exist")
})

# renameGear ----
test_that("renameGear works", {
    p <- example_manipulate_params
    # Get original gear names
    original_gears <- dimnames(p@selectivity)$gear
    gear1 <- original_gears[1]
    gear2 <- original_gears[2]

    # Define replacement
    replace <- c("new_gear1", "new_gear2")
    names(replace) <- c(gear1, gear2)

    # Rename gears
    p2 <- renameGear(p, replace)
    expect_false(identical(p2@time_modified, p@time_modified))

    # Check that gear_params is updated
    expect_true("new_gear1" %in% p2@gear_params$gear)
    expect_true("new_gear2" %in% p2@gear_params$gear)
    expect_false(gear1 %in% p2@gear_params$gear)
    expect_false(gear2 %in% p2@gear_params$gear)

    # Check that selectivity dimension names are updated
    new_gears <- dimnames(p2@selectivity)$gear
    expect_true("new_gear1" %in% new_gears)
    expect_true("new_gear2" %in% new_gears)
    expect_false(gear1 %in% new_gears)
    expect_false(gear2 %in% new_gears)

    # Check that catchability dimension names are updated
    expect_identical(dimnames(p2@catchability)$gear, new_gears)

    # Check that initial_effort names are updated
    expect_true("new_gear1" %in% names(p2@initial_effort))
    expect_true("new_gear2" %in% names(p2@initial_effort))
    expect_false(gear1 %in% names(p2@initial_effort))
    expect_false(gear2 %in% names(p2@initial_effort))

    # Check that the values in initial_effort are preserved
    expect_equal(p2@initial_effort[["new_gear1"]], p@initial_effort[[gear1]])
    expect_equal(p2@initial_effort[["new_gear2"]], p@initial_effort[[gear2]])

    # Check that params object is valid
    expect_true(validObject(p2))
})

test_that("renameGear warns on wrong names", {
    expect_error(renameGear(example_manipulate_params,
                            c(Trawler = "New_Trawl", NonExistent = "Other")),
                 "Trawler, NonExistent do not exist")
})

# expandSizeGrid ----
test_that("expandSizeGrid is deprecated", {
    expect_warning(expandSizeGrid(NS_params_small), "deprecated")
})

test_that("expandSizeGrid works", {
    params <- suppressWarnings(expandSizeGrid(NS_params_small, new_min_w = 1e-4, new_max_w = 50000))
    expect_lt(min(params@w), 1e-4)
    expect_gt(max(params@w), 50000)
    min_idx <- sum(params@w < min(NS_params_small@w)) + 1
    max_idx <- sum(params@w < max(NS_params_small@w)) + 1
    expect_equal(params@w[min_idx:max_idx], NS_params_small@w)
    expect_equal(params@search_vol[, min_idx:max_idx], NS_params_small@search_vol)
    expect_equal(params@metab[, min_idx:max_idx], NS_params_small@metab)
    expect_equal(params@initial_n[, min_idx:max_idx], NS_params_small@initial_n)
    expect_false(identical(params@time_modified, NS_params_small@time_modified))
})

test_that("expandSizeGrid preserves existing data", {
    params <- suppressWarnings(expandSizeGrid(NS_params_small, new_min_w = min(NS_params_small@w), new_max_w = max(NS_params_small@w)))
    expect_unchanged(params, NS_params_small)
})

test_that("expandSizeGrid preserves slots recently added to upgradeParams", {
    set.seed(42)
    params <- small_trait_grid_manipulate_params
    params <- randomise_upgrade_preserved_slots(params)

    p <- suppressWarnings(expandSizeGrid(params, new_max_w = max(params@w) * 10))
    old_sp <- seq_len(nrow(params@species_params))
    old_w <- seq_along(params@w)

    expect_equal(p@ext_diffusion[old_sp, old_w], params@ext_diffusion,
                 ignore_attr = TRUE)
    for (slot in c("initial_n_other", "other_dynamics", "other_encounter",
                   "other_mort", "other_params", "rates_funcs",
                   "use_predation_diffusion")) {
        expect_identical(slot(p, slot), slot(params, slot))
    }
})

test_that("expandSizeGrid preserves the class", {
    if (!methods::isClass("ExpandGridTestParams")) {
        methods::setClass("ExpandGridTestParams", contains = "MizerParams")
    }
    params <- as(small_trait_grid_manipulate_params, "ExpandGridTestParams")
    p <- suppressWarnings(expandSizeGrid(params, new_max_w = max(params@w) * 10))

    expect_s4_class(p, "ExpandGridTestParams")
    expect_true(validObject(p))
})

# adjustSizeGrid ----
test_that("adjustSizeGrid works for expansion and truncation", {
    # 1. No-op
    p_noop <- adjustSizeGrid(NS_params_small, new_min_w = min(NS_params_small@w), new_max_w = max(NS_params_small@w))
    expect_unchanged(p_noop, NS_params_small)

    # 2. Expansion (should match expandSizeGrid output)
    p_exp <- adjustSizeGrid(NS_params_small, new_min_w = 1e-4, new_max_w = 50000)
    p_exp_ref <- suppressWarnings(expandSizeGrid(NS_params_small, new_min_w = 1e-4, new_max_w = 50000))
    expect_equal(p_exp, p_exp_ref)

    # 3. Truncation
    # First, let's zero out abundance in the margins to avoid warnings
    params <- NS_params_small
    min_w_val <- params@w[5]
    max_w_val <- params@w[15]
    params@initial_n[, params@w < min_w_val] <- 0
    params@initial_n[, params@w > max_w_val] <- 0
    params@initial_n_pp[params@w_full < min_w_val] <- 0
    params@initial_n_pp[params@w_full > max_w_val] <- 0

    p_trunc <- expect_warning(adjustSizeGrid(params, new_min_w = min_w_val, new_max_w = max_w_val), NA)
    expect_equal(min(p_trunc@w), min_w_val)
    expect_equal(max(p_trunc@w), max_w_val)

    # 4. Truncation warning when abundance is present
    params <- NS_params_small
    # With original abundance (very high at small sizes), truncating at the bottom triggers species warning:
    expect_warning(adjustSizeGrid(params, new_min_w = params@w[3]),
                   "Non-negligible species biomass was lost due to grid truncation: Sprat")

    # Truncating at the top triggers species biomass warning:
    expect_warning(adjustSizeGrid(params, new_max_w = params@w[18]),
                   "Non-negligible species biomass was lost")

    # Truncating further down with low tol and non-zero resource at large sizes triggers resource biomass warning:
    params_pp <- params
    params_pp@resource_params$w_pp_cutoff <- 50000
    params_pp@initial_n_pp[] <- 1
    expect_warning(
        expect_warning(
            adjustSizeGrid(params_pp, new_max_w = params@w[18], tol = 1e-15),
            "Non-negligible species biomass was lost"
        ),
        "Non-negligible resource biomass"
    )

    # 5. Invalid parameters
    expect_error(adjustSizeGrid(NS_params_small, new_min_w = 10, new_max_w = 2),
                 "new_min_w must be smaller than new_max_w.")
    # Species w_max smaller than new_min_w
    expect_error(adjustSizeGrid(NS_params_small, new_min_w = 100),
                 "The following species have their maximum size w_max smaller than the new minimum size: Sprat")
    # Species w_min larger than new_max_w
    expect_error(adjustSizeGrid(NS_params_small, new_min_w = 0.00001, new_max_w = 0.0001),
                 "The following species have their minimum size w_min larger than the new maximum size:")
})

