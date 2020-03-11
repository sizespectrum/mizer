context("Backwards compatibility")

# The known values below were calculated with mizer version 1.0.1

# In version 1.0.1 it was still necessary to use data()
data(NS_species_params_gears)
data(inter)

test_that("MizerParams() works as in version 1", {
    params <- MizerParams(NS_species_params_gears, inter)
    expect_known_value(params@search_vol, "values/set_multispecies_model_search_vol")
    expect_known_value(params@intake_max, "values/set_multispecies_model_intake_max")
    expect_known_value(params@psi,        "values/set_multispecies_model_psi")
    #expect_known_value(params@std_metab + params@activity, "values/set_multispecies_model_metab")
    #expect_known_value(params@metab,      "values/set_multispecies_model_metab")
    expect_known_value(params@mu_b,      "values/set_multispecies_model_mu_b")
    expect_known_value(params@ft_pred_kernel_e, "values/set_multispecies_model_ft_pred_kernel_e")
    expect_known_value(params@ft_pred_kernel_p, "values/set_multispecies_model_ft_pred_kernel_p")
    expect_known_value(params@rr_pp, "values/set_multispecies_model_rr_pp")
    expect_known_value(params@cc_pp, "values/set_multispecies_model_cc_pp")
    expect_known_value(params@initial_n, "values/set_multispecies_model_initial_n")
    expect_known_value(params@initial_n_pp, "values/set_multispecies_model_initial_n_pp")
    expect_known_value(getPredRate(params, n = params@initial_n, n_pp = params@initial_n_pp), 
                       "values/set_multispecies_model_pred_rate")
    expect_known_value(getPhiPrey(params, n = params@initial_n, n_pp = params@initial_n_pp), 
                       "values/set_multispecies_model_phi_prey")
    sim <- project(params, t_max = 0.1, t_save = 0.1)
    expect_known_value(sim@n[2, , ],  "values/set_multispecies_model_n")
    expect_known_value(sim@n_pp[2, ], "values/set_multispecies_model_n_pp")
})

test_that("set_trait_model() works as in version 1", {
    params <- set_trait_model()
    expect_known_value(params@search_vol, "values/set_trait_model_search_vol")
    expect_known_value(params@intake_max, "values/set_trait_model_intake_max")
    expect_known_value(params@psi,        "values/set_trait_model_psi")
    #expect_known_value(params@std_metab + params@activity, "values/set_trait_model_metab")
    #expect_known_value(params@metab,      "values/set_trait_model_metab")
    expect_known_value(params@mu_b,      "values/set_trait_model_mu_b")
    expect_known_value(params@ft_pred_kernel_e, "values/set_trait_model_ft_pred_kernel_e")
    expect_known_value(params@ft_pred_kernel_p, "values/set_trait_model_ft_pred_kernel_p")
    expect_known_value(params@rr_pp, "values/set_trait_model_rr_pp")
    expect_known_value(params@cc_pp, "values/set_trait_model_cc_pp")
    expect_known_value(params@initial_n, "values/set_trait_model_initial_n")
    expect_known_value(params@initial_n_pp, "values/set_trait_model_initial_n_pp")
    expect_known_value(getPredRate(params, n = params@initial_n, n_pp = params@initial_n_pp), 
                       "values/set_trait_model_pred_rate")
    expect_known_value(getPhiPrey(params, n = params@initial_n, n_pp = params@initial_n_pp), 
                       "values/set_trait_phi_prey")
    sim <- project(params, t_max = 0.1, t_save = 0.1)
    expect_known_value(sim@n[2, , ],  "values/set_trait_model_n")
    expect_known_value(sim@n_pp[2, ], "values/set_trait_model_n_pp")
})

test_that("set_community_model() works as in version 1", {
    params <- set_community_model()
    expect_known_value(params@search_vol, "values/set_community_model_search_vol")
    expect_known_value(params@intake_max, "values/set_community_model_intake_max")
    expect_known_value(params@psi,        "values/set_community_model_psi")
    #expect_known_value(params@std_metab + params@activity, "values/set_community_model_metab")
    #expect_known_value(params@metab,      "values/set_community_model_metab")
    expect_known_value(params@mu_b,      "values/set_community_model_mu_b")
    expect_known_value(params@ft_pred_kernel_e, "values/set_community_model_ft_pred_kernel_e")
    expect_known_value(params@ft_pred_kernel_p, "values/set_community_model_ft_pred_kernel_p")
    expect_known_value(params@rr_pp, "values/set_community_model_rr_pp")
    expect_known_value(params@cc_pp, "values/set_community_model_cc_pp")
    expect_known_value(params@initial_n, "values/set_community_model_initial_n")
    expect_known_value(params@initial_n_pp, "values/set_community_model_initial_n_pp")
    expect_known_value(getPredRate(params, n = params@initial_n, n_pp = params@initial_n_pp), 
                       "values/set_community_model_pred_rate")
    expect_known_value(getPhiPrey(params, n = params@initial_n, n_pp = params@initial_n_pp), 
                       "values/set_community_phi_prey")
    sim <- project(params, t_max = 0.1, t_save = 0.1)
    expect_known_value(sim@n[2, , ],  "values/set_community_model_n")
    expect_known_value(sim@n_pp[2, ], "values/set_community_model_n_pp")
})
