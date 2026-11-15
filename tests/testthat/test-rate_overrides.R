test_that("override input is explicit and legacy whole-array freezes survive", {
    withr::local_options(mizer_info_level = 0)
    p <- NS_params_small
    original <- p
    metab(p) <- metab(p) * 0.5
    saved <- withr::local_tempfile(fileext = ".rds")
    saveRDS(p, saved)
    p <- readRDS(saved)
    expect_setequal(rate_overrides(p)$species, species_params(p)$species)
    p <- setMetabolicRate(p, species = "Cod", reset = TRUE)
    expect_equal(metab(p)["Cod", ], metab(original)["Cod", ], ignore_attr = TRUE)
    expect_setequal(rate_overrides(p)$species, c("Sprat", "Herring"))
    saveRDS(p, saved)
    expect_equal(rate_overrides(readRDS(saved)), rate_overrides(p))
    expect_error(setMetabolicRate(p, species = "unknown"), "unique species")
    expect_error(setMetabolicRate(p, species = c("Cod", "Cod")), "unique species")
    expect_error(setMetabolicRate(p, species = character()), "unique species")
    expect_error(setMetabolicRate(p, metab = 1, species = "Cod"), "Supply an array")
    expect_error(setMetabolicRate(p, metab = rep(NA_real_, length(w(p))),
                                  species = "Cod"), "finite")
    rows <- metab(p)[c("Cod", "Sprat"), , drop = FALSE]
    expect_error(setMetabolicRate(p, metab = rows,
                                  species = c("Sprat", "Cod")), "supplied order")
    q <- setMetabolicRate(original, metab = rows, species = c("Cod", "Sprat"))
    expect_setequal(rate_overrides(q)$species, c("Cod", "Sprat"))
    metab(q) <- metab(q)
    expect_equal(nrow(rate_overrides(q)), 3L)
    q <- `metab<-`(original, metab(original))
    expect_equal(nrow(rate_overrides(q)), 3L)
    q <- `catchability<-`(original, catchability(original))
    expect_equal(nrow(rate_overrides(q)), 3L)
})

test_that("override reporting retains inputs and identifies affected rates", {
    p <- NS_params_small
    before <- p@given_species_params
    metab(p, species = "Cod") <- rep(0.2, length(w(p)))
    report <- rate_overrides(p)
    expect_equal(report$rate, "metab")
    expect_equal(report$species, "Cod")
    expect_true("ks" %in% report$parameters[[1]])
    expect_match(report$reset, 'species = "Cod"', fixed = TRUE)
    expect_equal(p@given_species_params, before)
    expect_equal(p@species_params, NS_params_small@species_params)
    expect_output(print(species_params(p)), "metab for Cod")
    expect_output(print(species_params(p)[, "ks", drop = FALSE]), "metab for Cod")
    expect_false(any(grepl("metab for Cod", capture.output(
        print(species_params(p)["Herring", ])))))
    # Evaluate the documented reset call and recover the original recipe.
    params <- p
    params <- eval(parse(text = report$reset))
    expect_equal(metab(params), metab(NS_params_small), ignore_attr = TRUE)
})
