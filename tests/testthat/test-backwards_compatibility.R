# The known values below were calculated with mizer version 1.0.1

# In version 1.0.1 it was still necessary to use data()
data(NS_species_params_gears)
data(inter)

test_that("MizerParams() works as in version 2.5.1", {
  old <- getOption("mizer_defaults_edition")
  on.exit(options(mizer_defaults_edition = old), add = TRUE)
  options(mizer_defaults_edition = 1)
  # expect_warning(params <- MizerParams(NS_species_params_gears, inter), "deprecated")
  # warning no longer thrown - NS_species_params_gears is now 2.5.1 (see `params@mizer_version`)
  params <- newMultispeciesParams(NS_species_params_gears, inter, info_level = 0)
  # expect_known_value(params@search_vol, "values/set_multispecies_model_search_vol")
  # expect_known_value(params@intake_max, "values/set_multispecies_model_intake_max")
  # expect_known_value(params@psi,        "values/set_multispecies_model_psi")
  #expect_known_value(params@std_metab + params@activity, "values/set_multispecies_model_metab")
  # expect_known_value(params@metab,      "values/set_multispecies_model_metab")
  # expect_known_value(params@mu_b,      "values/set_multispecies_model_mu_b")
  # expect_known_value(params@rr_pp, "values/set_multispecies_model_rr_pp")
  # expect_known_value(params@cc_pp, "values/set_multispecies_model_cc_pp")
  # expect_known_value(params@initial_n, "values/set_multispecies_model_initial_n")
  # expect_known_value(params@initial_n_pp, "values/set_multispecies_model_initial_n_pp")

  expect_snapshot(params@search_vol)
  expect_snapshot(params@intake_max)
  expect_snapshot(params@psi)
  expect_snapshot(params@metab)
  expect_snapshot(params@mu_b)
  expect_snapshot(params@rr_pp)
  expect_snapshot(params@cc_pp)
  expect_snapshot(params@initial_n)
  expect_snapshot(params@initial_n_pp)

  # The predation rate is different because in the old version the predation
  # kernel was cut off at beta + 3 sigma
  # new <- getPredRate(params)
  # df <- melt(new) %>% mutate(type = "new")
  # old <- readRDS("values/set_multispecies_model_pred_rate")
  # dimnames(old) <- dimnames(new)
  # dfo <- melt(old) %>% mutate(type = "old")
  # dfs <- rbind(df, dfo)
  # dfs$sp <- as.factor(dfs$sp)
  # ggplot(dfs) +
  #     geom_line((aes(x = w_prey, y = value, colour = sp, linetype = type))) +
  #     scale_x_log10() +
  #     scale_y_log10()
  # max(abs(old - new))
  # expect_warning(gpp <- getPhiPrey(params, n = initialN(params),
  #                                  n_pp = initialNResource(params)),
  #                "deprecated")
  # expect_known_value(gpp, "values/set_multispecies_model_phi_prey",
  #                    check.attributes = FALSE,
  #                    tolerance = 0.0005)
  expect_snapshot(getPhiPrey(params, n = initialN(params), n_pp = initialNResource(params)))
  # The discrepancy is small, see the following graph:
  # new <- getPhiPrey(params, n = initialN(params), n_pp = initialNResource(params))
  # df <- melt(new) %>% mutate(type = "new")
  # old <- readRDS("values/set_multispecies_model_phi_prey")
  # dimnames(old) <- dimnames(new)
  # dfo <- melt(old) %>% mutate(type = "old")
  # dfs <- rbind(df, dfo)
  # dfs$sp <- as.factor(dfs$sp)
  # ggplot(dfs) +
  #     geom_line((aes(x = w, y = value, colour = sp, linetype = type))) +
  #     scale_x_log10() +
  #     scale_y_log10()
  # max(abs(old - new))
})

test_that("set_trait_model() works as in version 1", {
  old <- getOption("mizer_defaults_edition")
  on.exit(options(mizer_defaults_edition = old), add = TRUE)
  options(mizer_defaults_edition = 1)
  # expect_warning(params <- set_trait_model(), "deprecated")
  # warning no longer thrown - default params are now 2.5.1 (see `params@mizer_version`)
  params <- newTraitParams()
  # expect_known_value(params@search_vol, "values/set_trait_model_search_vol")
  # expect_known_value(params@intake_max, "values/set_trait_model_intake_max")
  # expect_known_value(params@psi,        "values/set_trait_model_psi")
  #expect_known_value(params@std_metab + params@activity, "values/set_trait_model_metab")
  # expect_known_value(params@metab,      "values/set_trait_model_metab")
  # expect_known_value(params@mu_b,      "values/set_trait_model_mu_b")
  # expect_known_value(params@rr_pp, "values/set_trait_model_rr_pp")
  # expect_known_value(params@cc_pp, "values/set_trait_model_cc_pp")
  # expect_known_value(params@initial_n, "values/set_trait_model_initial_n")
  # expect_known_value(params@initial_n_pp, "values/set_trait_model_initial_n_pp")

  expect_snapshot(params@search_vol)
  expect_snapshot(params@intake_max)
  expect_snapshot(params@psi)
  expect_snapshot(params@metab)
  expect_snapshot(params@mu_b)
  expect_snapshot(params@rr_pp)
  expect_snapshot(params@cc_pp)
  expect_snapshot(params@initial_n)
  expect_snapshot(params@initial_n_pp)

  # The predation rate is different because in the old version the predation
  # kernel was cut off at beta + 3 sigma
  # new <- getPredRate(params)
  # df <- melt(new) %>% mutate(type = "new")
  # old <- readRDS("values/set_trait_model_pred_rate")
  # dimnames(old) <- dimnames(new)
  # dfo <- melt(old) %>% mutate(type = "old")
  # dfs <- rbind(df, dfo)
  # dfs$sp <- as.factor(dfs$sp)
  # ggplot(dfs) +
  #     geom_line((aes(x = w_prey, y = value, colour = sp, linetype = type))) +
  #     scale_x_log10() +
  #     scale_y_log10()
  # expect_known_value(getPhiPrey(params, n = initialN(params), n_pp = initialNResource(params)),
  #                    "values/set_trait_phi_prey",
  #                    check.attributes = FALSE,
  #                    tolerance = 1e-5) # set_trait_phi_prey doesn't have dimnames, needs updating
  expect_snapshot(getPhiPrey(params, n = initialN(params), n_pp = initialNResource(params)))

  # The discrepancy is too small to see in the following graph:
  # new <- getPhiPrey(params, n = initialN(params), n_pp = initialNResource(params))
  # df <- melt(new) %>% mutate(type = "new")
  # old <- readRDS("values/set_trait_phi_prey")
  # dimnames(old) <- dimnames(new)
  # dfo <- melt(old) %>% mutate(type = "old")
  # dfs <- rbind(df, dfo)
  # dfs$sp <- as.factor(dfs$sp)
  # ggplot(dfs) +
  #     geom_line((aes(x = w, y = value, colour = sp, linetype = type))) +
  #     scale_x_log10() +
  #     scale_y_log10()
  # max(abs(old - new))
})

test_that("set_trait_model uses documented grid arguments", {
  params <- set_trait_model(max_w = 2e5, min_w_pp = 1e-8) |>
      suppressMessages() |> suppressWarnings()
  expect_equal(max(w(params)), 2e5)
  expect_equal(min(w_full(params)), 1e-8)
})

test_that("set_community_model() works as in version 1", {
  old <- getOption("mizer_defaults_edition")
  on.exit(options(mizer_defaults_edition = old), add = TRUE)
  options(mizer_defaults_edition = 1)
  # expect_warning(params <- set_community_model(), "deprecated")
  # warning no longer thrown - default params are now 2.5.1 (see `params@mizer_version`)
  params <- newCommunityParams()
  # expect_known_value(params@search_vol, "values/set_community_model_search_vol")
  # expect_known_value(params@intake_max, "values/set_community_model_intake_max")
  # expect_known_value(params@psi,        "values/set_community_model_psi")
  #expect_known_value(params@std_metab + params@activity, "values/set_community_model_metab")
  # expect_known_value(params@metab,      "values/set_community_model_metab")
  # expect_known_value(params@mu_b,      "values/set_community_model_mu_b")
  # expect_known_value(params@rr_pp, "values/set_community_model_rr_pp")
  # expect_known_value(params@cc_pp, "values/set_community_model_cc_pp")
  # expect_known_value(params@initial_n, "values/set_community_model_initial_n")
  # expect_known_value(params@initial_n_pp, "values/set_community_model_initial_n_pp")

  expect_snapshot(params@search_vol)
  expect_snapshot(params@intake_max)
  expect_snapshot(params@psi)
  expect_snapshot(params@metab)
  expect_snapshot(params@mu_b)
  expect_snapshot(params@rr_pp)
  expect_snapshot(params@cc_pp)
  expect_snapshot(params@initial_n)
  expect_snapshot(params@initial_n_pp)

  # The predation rate is different because in the old version the predation
  # kernel was cut off at beta + 3 sigma, see following graph
  # new <- getPredRate(params)
  # df <- melt(new) %>% mutate(type = "new")
  # old <- readRDS("values/set_community_model_pred_rate")
  # dimnames(old) <- dimnames(new)
  # dfo <- melt(old) %>% mutate(type = "old")
  # dfs <- rbind(df, dfo)
  # dfs$sp <- as.factor(dfs$sp)
  # ggplot(dfs) +
  #     geom_line((aes(x = w_prey, y = value, colour = sp, linetype = type))) +
  #     scale_x_log10() +
  #     scale_y_log10()
  # max(abs(new - old))
  # expect_known_value(getPhiPrey(params, n = initialN(params), n_pp = initialNResource(params)),
  #                    "values/set_community_phi_prey",
  #                    check.attributes = FALSE,
  #                    tolerance = 42)
  expect_snapshot(getPhiPrey(params, n = initialN(params), n_pp = initialNResource(params)))
  # The discrepancy is small, see the following graph:
  # new <- getPhiPrey(params, n = initialN(params), n_pp = initialNResource(params))
  # df <- melt(new) %>% mutate(type = "new")
  # old <- readRDS("values/set_community_phi_prey")
  # dimnames(old) <- dimnames(new)
  # dfo <- melt(old) %>% mutate(type = "old")
  # dfs <- rbind(df, dfo)
  # dfs$sp <- as.factor(dfs$sp)
  # ggplot(dfs) +
  #     geom_line((aes(x = w, y = value, colour = sp, linetype = type))) +
  #     scale_x_log10() +
  #     scale_y_log10()
  # max(abs(new - old))
})

test_that("set_community_model sets constant-community behaviour", {
  params <- set_community_model() |>
      suppressMessages() |> suppressWarnings()
  expect_identical(params@rates_funcs$RDD, "constantRDD")
  expect_true(all(params@psi == 0))
  expect_true(all(is.na(params@species_params$w_mat)))
})

test_that("set_multispecies_model keeps documented legacy defaults", {
  expect_warning(
    suppressMessages(params <- set_multispecies_model(NS_species_params_gears, inter, info_level = 0)),
    "deprecated"
  )
  expect_equal(max(w(params)), 1.1 * max(NS_species_params_gears$w_max))

  reduced <- NS_species_params_gears[
    ,
    setdiff(
      names(NS_species_params_gears),
      c("gear", "alpha", "erepro", "sel_func", "knife_edge_size",
        "catchability", "gamma", "ks", "m")
    )
  ]
  expect_warning(
    suppressMessages(params2 <- set_multispecies_model(reduced, inter, info_level = 0)),
    "deprecated"
  )
  sp <- species_params(params2)
  expect_identical(sp$gear, sp$species)
  expect_true(all(sp$alpha == 0.6))
  expect_true(all(sp$erepro == 1))
  expect_true(all(sp$sel_func == "knife_edge"))
  expect_identical(sp$knife_edge_size, sp$w_mat)
  expect_true(all(sp$catchability == 1))
  expect_true(all(sp$m == 1))
  expect_equal(sp$ks, sp$h * 0.2)
})

test_that("getPhiPrey deprecates softly and matches encounter over search volume", {
  params <- NS_params_small
  expect_warning(
    phi <- getPhiPrey(params, n = initialN(params), n_pp = initialNResource(params)),
    "deprecated"
  )
  expect_equal(phi, getEncounter(params) / getSearchVolume(params))
})
