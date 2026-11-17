# Shared plotting fixtures; unrelated tests do not force these projections.
# Numeric species names and the 30-bin grid preserve existing visual baselines.
delayedAssign("plot_params_small", with_defaults_edition(1, {
    sp <- NS_species_params_gears_small
    sp$species <- seq_len(nrow(sp))
    sp$pred_kernel_type <- "truncated_lognormal"
    newMultispeciesParams(sp, inter_small, no_w = 30,
                          n = 2 / 3, p = 0.7, lambda = 2.8 - 2 / 3,
                          info_level = 0)
}))
delayedAssign("plot_sim_small", project(plot_params_small, effort = 1,
                                        t_max = 3, dt = 1, t_save = 1,
                                        progress_bar = FALSE))
delayedAssign("plot_sim0_small", project(plot_params_small, effort = 0,
                                         t_max = 3, dt = 1, t_save = 1,
                                         progress_bar = FALSE))
delayedAssign("plot_params_background_small",
              markBackground(plot_params_small, species = 1))
delayedAssign("plot_params_single_small", with_defaults_edition(1, {
    newMultispeciesParams(data.frame(species = 1, w_max = 1000, h = 30),
                          no_w = 30, info_level = 0)
}))

expect_doppelganger <- function(title, fig, ...) {
    testthat::skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger(title, fig, ...)
}

expect_ggplot <- function(fig) {
    testthat::expect_s3_class(fig, "ggplot")
}

tooltip_fields <- function(gp, trace = 1) {
    tip <- gp$x$data[[trace]]$text[[1]]
    trimws(sub(":.*", "", strsplit(tip, "<br />", fixed = TRUE)[[1]]))
}
