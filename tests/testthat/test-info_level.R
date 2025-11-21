test_that("info_level works in newCommunityParams", {
    expect_message(newCommunityParams(info_level = 3), "Using f0, h, lambda, kappa")
    expect_silent(newCommunityParams(info_level = 0))
})

test_that("info_level works in newTraitParams", {
    expect_message(newTraitParams(info_level = 3), "Using f0, h, lambda, kappa")
    expect_silent(newTraitParams(info_level = 0))
})

test_that("info_level works in newMultispeciesParams", {
    sp <- NS_species_params
    expect_message(newMultispeciesParams(sp, info_level = 3), "Because")
    expect_silent(newMultispeciesParams(sp, info_level = 0))
})

test_that("info_level works in setParams", {
    params <- NS_params
    # setParams is called by newMultispeciesParams, but we can call it directly on a params object
    # However, setParams usually takes a params object that is already valid or being built.
    # Let's try modifying a param that triggers a default calculation.
    # But setParams doesn't recalculate defaults for existing slots unless they are missing?
    # Actually setParams calls setMaxIntakeRate etc.
    # If we set h to NA, it might trigger a message.
    params@species_params$h <- NA
    expect_message(setParams(params, info_level = 3), "Because")
    expect_silent(setParams(params, info_level = 0))
})

test_that("info_level works in projectToSteady", {
    params <- newTraitParams(no_sp = 2, info_level = 0)
    expect_message(projectToSteady(params, t_per = 1, t_max = 1, info_level = 3), "did not converge")
    expect_silent(projectToSteady(params, t_per = 1, t_max = 1, info_level = 0))
})

test_that("info_level works in matchYields", {
    params <- NS_params
    species_params(params)$yield_observed <- NA
    # matchYields is deprecated, suppress warning
    suppressWarnings({
        expect_message(matchYields(params, info_level = 3), "The following species have no yield observations")
        expect_silent(matchYields(params, info_level = 0))
    })
})
