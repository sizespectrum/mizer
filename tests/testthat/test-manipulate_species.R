# addSpecies ----
test_that("addSpecies works when adding a second identical species", {
    p <- newTraitParams()
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
    p <- NS_params
    species_params <- p@species_params[5, ]
    expect_error(addSpecies(p, species_params),
                 "You can not add species that are already there.")
})
test_that("addSpecies handles gear params correctly", {
    p <- newTraitParams(no_sp = 2)
    sp <- data.frame(species = c("new1", "new2"),
                     w_max = c(10, 100),
                     k_vb = c(4, 1),
                     n = 2 / 3,
                     p = 2 / 3)
    gp <- data.frame(gear = c("gear1", "gear2", "gear1"),
                     species = c("new1", "new2", "new2"),
                     sel_func = "knife_edge",
                     knife_edge_size = c(5, 5, 50))

    # If no initial_effort for new gear is provided, it is 0
    # Wrapping in `expect_warning()` to ignore warnings about unrealistic
    # reproductive efficiency
    (pa <- addSpecies(p, sp, gp)) |>
        expect_message() |>
        expect_warning()
    expect_identical(pa@initial_effort,
                     c(knife_edge_gear = 0, gear1 = 0, gear2 = 0))
    expect_identical(nrow(pa@gear_params), 5L)

    # effort for existing gear is not changed
    extra_effort <- c(gear1 = 2, gear2 = 3)
    (pa <- addSpecies(p, sp, gp, initial_effort = extra_effort)) |>
        expect_message() |>
        expect_warning()
    expect_identical(pa@initial_effort, c(knife_edge_gear = 0, extra_effort))

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
    p <- newTraitParams(no_sp = 2)
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
    params <- NS_params
    # change a slot to test that such changes will be preserved
    params <- setMaxIntakeRate(params, 2 * getMaxIntakeRate(params))

    (p <- addSpecies(params, sp)) |>
        expect_message()
    expect_identical(p@w[1:100], params@w)
    expect_identical(p@w_full[seq_along(params@w_full)], params@w_full)
    expect_lte(5e4, max(p@w))
    # changed rates are preserved
    expect_equal(getMaxIntakeRate(p)[1:12, 1:100],
                 getMaxIntakeRate(params), ignore_attr = TRUE)
})
test_that("addSpecies works when adding a species with a smaller w_min", {
    sp <- data.frame(species = "Blue whale", w_max = 5e4, w_min = 1e-5,
                     w_mat = 1e3, beta = 1000, sigma = 2,
                     k_vb = 0.6, gear = 'Whale hunter')
    params <- NS_params
    # change a slot to test that such changes will be preserved
    params <- setMaxIntakeRate(params, 2 * getMaxIntakeRate(params))

    (p <- addSpecies(params, sp)) |>
        expect_message()
    expect_equal(p@w[28:127], params@w)
    expect_equal(p@w_full[seq_along(params@w_full)], params@w_full)
    expect_gte(1e-5, min(p@w))
    # changed rates are preserved
    expect_equal(getMaxIntakeRate(p)[1:12, 28:127],
                 getMaxIntakeRate(params), ignore_attr = TRUE)
})

test_that("addSpecies has other documented properties", {
    sp <- data.frame(species = c("new1", "new2"),
                     w_max = c(10, 100),
                     k_vb = c(4, 1),
                     n = 2 / 3,
                     p = 2 / 3)
    (p <- addSpecies(NS_params, sp)) |>
        expect_message()

    # New species have 0 reproduction level
    expect_equal(getReproductionLevel(p)[13:14],
                 c(new1 = 1 / 4, new2 = 1 / 4))

    # Maximum of ratio between new species density and Sheldon density is 1/100
    fraction <- p@initial_n[13, ] /
        (p@resource_params$kappa * p@w ^ -p@resource_params$lambda)
    expect_equal(max(fraction), 1 / 100)
})

test_that("Added species stay at low abundance", {
    # Use example from man page
    params <- newTraitParams()
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

    params <- newTraitParams()
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

# removeSpecies ----
test_that("removeSpecies works", {
    remove <- NS_species_params$species[2:11]
    reduced <- NS_species_params[!(NS_species_params$species %in% remove), ]
    params <- newMultispeciesParams(NS_species_params, no_w = 20,
                                    max_w = 39900, min_w_pp = 9e-14,
                                    info_level = 0)
    p1 <- removeSpecies(params, species = remove)
    expect_equal(nrow(p1@species_params), nrow(params@species_params) - 10)
    p2 <- newMultispeciesParams(reduced, no_w = 20,
                                max_w = 39900, min_w_pp = 9e-14, info_level = 0)
    p2@linecolour[2] = "#a08dfb" # update line colour
    expect_equal(p1, p2, ignore_attr = TRUE)
    sim1 <- project(p1, t_max = 0.4, t_save = 0.4)
    sim2 <- project(p2, t_max = 0.4, t_save = 0.4)
    expect_identical(sim1@n[2, 2, ], sim2@n[2, 2, ])
})
test_that("removeSpecies works with 3d pred kernel", {
    # It should make no difference whether we first set full pred kernel and
    # then remove a species, or the other way around.
    params1 <- NS_params
    params1 <- setPredKernel(params1, pred_kernel = getPredKernel(params1))
    params1 <- removeSpecies(params1, "Cod")
    params2 <- NS_params
    params2 <- removeSpecies(params2, "Cod")
    params2 <- setPredKernel(params2, pred_kernel = getPredKernel(params2))
    expect_unchanged(params1, params2)
})
test_that("removeSpecies works correctly on gear_params", {
    # We'll check that the resulting gear_params lead to the same selectivity
    # and catchability
    params <- removeSpecies(NS_params, "Cod")
    expect_equal(nrow(params@gear_params), 11)
    params2 <- setFishing(params)
    expect_unchanged(params, params2)
})

test_that("adding and then removing species leaves params unaltered", {
    params <- newMultispeciesParams(NS_species_params, info_level = 0)
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
    # comment on w_min_idx are not preserved
    comment(params@w_min_idx) <- NULL
    expect_unchanged(params, params2)
})

# renameSpecies ----
test_that("renameSpecies works", {
    sp <- NS_species_params
    p <- newMultispeciesParams(sp, info_level = 0)
    sp$species <- tolower(sp$species)
    replace <- NS_species_params$species
    names(replace) <- sp$species
    p2 <- newMultispeciesParams(sp, info_level = 0)
    p2 <- renameSpecies(p2, replace)
    p2@time_modified <- p@time_modified
    p2@time_created <- p@time_created
    expect_identical(p, p2)
})
test_that("renameSpecies warns on wrong names", {
    expect_error(renameSpecies(NS_params, c(Kod = "cod", Hadok = "haddock")),
                 "Kod, Hadok do not exist")
})

# expandSizeGrid ----
test_that("expandSizeGrid works", {
    params <- expandSizeGrid(NS_params, new_min_w = 1e-4, new_max_w = 50000)
    expect_lt(min(params@w), 1e-4)
    expect_gt(max(params@w), 50000)
})
