params <- NS_params
sim <- project(params, t_max = 0.1, t_save = 0.1)

## upgradeParams ----
test_that("upgradeParams leaves new params unchanged", {
    expect_equal(upgradeParams(params), params)
})
test_that("upgradeParams preserves comments", {
    comment(params) <- "test"
    for (slot in (slotNames(params))) {
        comment(slot(params, slot)) <- slot
    }
    expect_equal(upgradeParams(params), params)
})

## upgradeSim ----
test_that("upgradeSim leaves new sim unchanged", {
    expect_equal(upgradeSim(sim), sim)
})
test_that("upgradeSim preserves comments", {
    comment(sim) <- "test"
    for (slot in (slotNames(sim))) {
        comment(slot(sim, slot)) <- slot
    }
    expect_equal(upgradeSim(sim), sim)
})

test_that("Object from version 0.4 can be upgraded", {
    simc.0.4 <- readRDS("assets/simc.0.4.rds")
    sim <- upgradeSim(simc.0.4)
    expect_true(validObject(sim))
})
test_that("Object from version 1.0 can be upgraded", {
    simc.1.0 <- readRDS("assets/simc.1.0.rds")
    sim <- upgradeSim(simc.1.0)
    expect_true(validObject(sim))
})
