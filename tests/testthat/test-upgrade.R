params <- NS_params
sim <- project(params, t_max = 0.1, t_save = 0.1)

## upgradeParams ----
test_that("upgradeParams leaves new params unchanged", {
    expect_identical(upgradeParams(params), params)
})
test_that("upgradeParams preserves comments", {
    comment(params) <- "test"
    for (slot in (slotNames(params))) {
        comment(slot(params, slot)) <- slot
    }
    expect_identical(upgradeParams(params), params)
})

## upgradeSim ----
test_that("upgradeSim leaves new sim unchanged", {
    expect_identical(upgradeSim(sim), sim)
})
test_that("upgradeSim preserves comments", {
    comment(sim) <- "test"
    for (slot in (slotNames(sim))) {
        comment(slot(sim, slot)) <- slot
    }
    expect_identical(upgradeSim(sim), sim)
})

test_that("Object from version 0.4 can be upgraded", {
    simc.0.4 <- readRDS("assets/simc.0.4.rds")
    (sim <- validSim(simc.0.4)) |>
        expect_message("Initial effort has been set to 0") |>
        expect_warning("Your MizerSim object was created with an earlier")
    
    expect_true(validObject(sim))
})
test_that("Object from version 1.0 can be upgraded", {
    simc.1.0 <- readRDS("assets/simc.1.0.rds")
    (sim <- validSim(simc.1.0)) |>
        expect_message("Initial effort has been set to 0") |>
        expect_warning("Your MizerSim object was created with an earlier")
    expect_true(validObject(sim))
})
test_that("r_max is renamed", {
    params <- NS_params
    params@species_params$r_max <- params@species_params$R_max
    params@mizer_version <- "1.1"
    expect_message(upgradeParams(params),
                   "The 'r_max' column has been renamed to 'R_max'.")
})

test_that("Some functions work with params from earlier versions", {
    params.0.4 <- readRDS("assets/simc.0.4.rds")@params
    expect_warning(getEGrowth(params.0.4) |> expect_message(),
                   "Your MizerParams object was created with an earlier")
    expect_warning(plotFeedingLevel(params.0.4) |> expect_message(),
                   "Your MizerParams object was created with an earlier")
    expect_warning(project(params.0.4, t_max = 0.1) |> expect_message(),
                   "Your MizerParams object was created with an earlier")
    # renaming of resource dynamics functions
    slot(params.0.4, "srr", check = FALSE) <- "srrNone"
    (p4 <- suppressWarnings(upgradeParams(params.0.4))) |> expect_message()
    expect_identical(p4@rates_funcs$RDD, "noRDD")
})
