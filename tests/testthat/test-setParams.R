context("Setting parameters")
data(NS_species_params_gears)
data(NS_species_params)
data(inter)
no_sp <- nrow(NS_species_params)


## setResourceEncounter ----
test_that("setResourceEncounter works", {
    resource_dynamics <-
        list("detritus" = function(params, n, n_pp, B, rates, dt, ...) B["detritus"],
             "carrion" = function(params, n, n_pp, B, rates, dt, ...) B["carrion"])
    species_params <- NS_species_params
    species_params$rho_detritus <- 1:no_sp
    species_params$rho_carrion <- no_sp:1
    params <- MizerParams(species_params, resource_dynamics = resource_dynamics)
    expect_equal(params@rho[2, 1, 1], 2 * params@w[1]^params@n)
    expect_equal(params@rho[2, 2, 1], (no_sp - 1) * params@w[1]^params@n)
})