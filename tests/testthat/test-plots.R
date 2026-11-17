# Plot fixtures use edition 1, as do the existing visual baselines.
withr::local_options(mizer_defaults_edition = 1)
params <- plot_params_small
sim <- plot_sim_small
species <- c(1, 2)

test_that("the summary plots route the spectrum arguments to plotSpectra only", {
    # `per_log_size` is an error on a rate that is not a density, so passing it
    # to every panel made `plot(sim, per_log_size = TRUE)` fail even though the
    # equivalent `power = 2` worked.
    expect_identical(names(shared_plot_args(list(species = species, power = 2,
                                                 biomass = TRUE,
                                                 per_log_size = TRUE,
                                                 resource = FALSE))),
                     "species")
    # Arguments the other panels do understand are still shared
    expect_identical(names(shared_plot_args(list(total = TRUE,
                                                 background = FALSE))),
                     c("total", "background"))
    # Unnamed dots are not dropped
    expect_length(shared_plot_args(list(1, 2)), 2)
    expect_length(shared_plot_args(list()), 0)

    withr::local_pdf(NULL)
    expect_silent(plot(sim, species = species, per_log_size = TRUE))
    expect_silent(plot(sim, species = species, biomass = TRUE,
                       per_log_size = TRUE))
    expect_silent(plot(sim, species = species, power = 2))
    expect_silent(plot(params, species = species, per_log_size = TRUE))
})
