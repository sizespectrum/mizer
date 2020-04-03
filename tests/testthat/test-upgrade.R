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
# TODO: add tests that upgrade objects saved from previous versions.

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
