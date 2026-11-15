# Growth curves and their interactive compatibility wrapper.

#' Plot growth curves
#'
#' `r lifecycle::badge("experimental")`
#' The growth curves represent the average age of all the living fish of a
#' species as a function of their size. So it would be natural to plot size
#' on the x-axis. But to follow the usual convention from age-based models, we
#' plot size on the y-axis and age on the x-axis.
#'
#' In each panel for a single species, a horizontal line is included that
#' indicate the maturity size of the species and a vertical line indicating its
#' maturity age.
#'
#' If size at age data is passed via the `size_at_age` argument, this is plotted
#' on top of the growth curve. When comparing this to the growth curves, you
#' need to remember that the growth curves should only represent the average
#' age at each size. So a scatter in the x-direction around the curve is to be
#' expected.
#'
#' If the species parameters contain the variables `a` and `b` for length to
#' weight conversion and the von Bertalanffy parameter `k_vb`, `w_inf` (and
#' optionally `t0`), then the von Bertalanffy growth curve is superimposed in
#' black. Note that the von Bertalanffy curve (which approximates the average
#' length at each age) should not be compared directly to the mizer growth
#' curves (which approximate the average age at each length).
#'
#' @param object An object of class \linkS4class{MizerSim} or
#'   \linkS4class{MizerParams}.
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param max_age The age up to which to run the growth curve. Default is 20.
#' @param percentage Boolean value. If TRUE, the size is given as a percentage
#'   of the maximal size.
#' @param species_panel If TRUE (default), and `percentage = FALSE`, display all
#'   species as facets. Otherwise puts all species into a single panel.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param size_at_age A data frame with observed size at age data to be plotted
#'   on top of growth curve graphs. Should contain columns `species` (species
#'   name as used in the model), `age` (in years) and either `weight` (in grams)
#'   or `length` (in cm). If both `weight` and `length` are provided, only
#'   `weight` is used.
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default is FALSE.
#' @param log_x If `TRUE`, use a log10 x-axis. Default is `FALSE`.
#' @param log_y If `TRUE`, use a log10 y-axis. Default is `FALSE`.
#' @param log Character string specifying which axes should use log10 scales,
#'   in the same form as the base [plot()] argument. For example, `"x"`,
#'   `"y"`, `"xy"` or `""`. If supplied, this overrides `log_x` and `log_y`.
#' @param ... Unused.
#' @return A ggplot2 object
#' @export
#' @family plotting functions
#' @seealso [plotting_functions]
#' @examples
#' \donttest{
#' params <-  NS_params
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotGrowthCurves(sim, percentage = TRUE)
#' plotGrowthCurves(sim, species = "Cod", max_age = 24)
#' plotGrowthCurves(sim, species_panel = TRUE)
#'
#' # Returning the data frame
#' fr <- plotGrowthCurves(sim, return_data = TRUE)
#' str(fr)
#' }
#' @rdname plotGrowthCurves
#' @export
plotGrowthCurves <- function(object, species = NULL, max_age = 20,
                             percentage = FALSE, species_panel = FALSE,
                             highlight = NULL, size_at_age = NULL,
                             return_data = FALSE, log_x = FALSE,
                             log_y = FALSE, log = NULL, ...) {
    UseMethod("plotGrowthCurves")
}

#' @rdname plotGrowthCurves
#' @usage NULL
#' @export
plotGrowthCurves.MizerSim <- function(object, species = NULL,
                             max_age = 20, percentage = FALSE,
                             species_panel = FALSE, highlight = NULL,
                             size_at_age = NULL,
                             return_data = FALSE,
                             log_x = FALSE, log_y = FALSE, log = NULL, ...) {
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    assert_that(is.flag(percentage),
                is.flag(species_panel),
                is.flag(return_data),
                is.number(max_age))
    params <- object@params
    params <- finalParams(object)
    plot_growth_curves(params, species = species,
                       max_age = max_age, percentage = percentage,
                       species_panel = species_panel, highlight = highlight,
                       log_x = log_axes$log_x, log_y = log_axes$log_y,
                       size_at_age = size_at_age,
                       return_data = return_data)
}

#' @rdname plotGrowthCurves
#' @usage NULL
#' @export
plotGrowthCurves.MizerParams <- function(object, species = NULL,
                             max_age = 20, percentage = FALSE,
                             species_panel = FALSE, highlight = NULL,
                             size_at_age = NULL,
                             return_data = FALSE,
                             log_x = FALSE, log_y = FALSE, log = NULL, ...) {
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    assert_that(is.flag(percentage),
                is.flag(species_panel),
                is.flag(return_data),
                is.number(max_age))
    params <- validParams(object)
    plot_growth_curves(params, species = species,
                       max_age = max_age, percentage = percentage,
                       species_panel = species_panel, highlight = highlight,
                       log_x = log_axes$log_x, log_y = log_axes$log_y,
                       size_at_age = size_at_age,
                       return_data = return_data)
}

#' Build the growth-curves plot
#'
#' Internal worker shared by the [plotGrowthCurves()] methods. It computes the
#' modelled size at age, optionally adds a von Bertalanffy curve and observed
#' size-at-age data, and draws the plot.
#'
#' @param params A MizerParams object.
#' @param species The species to be plotted.
#' @param max_age The age up to which to plot the growth curve.
#' @param percentage If `TRUE`, size is shown as a percentage of maximum size.
#' @param species_panel If `TRUE`, each species is shown in its own panel.
#' @param highlight Name or vector of names of species to be highlighted.
#' @param log_x,log_y Logical flags for log10 axes.
#' @param size_at_age Optional data frame of observed size-at-age data.
#' @param return_data If `TRUE`, return the plotting data frame instead of the
#'   plot.
#' @return A `mizer_plot` (ggplot2) object, or the plotting data frame if
#'   `return_data = TRUE`.
#' @keywords internal
plot_growth_curves <- function(params, species,
                               max_age, percentage,
                               species_panel, highlight,
                               log_x, log_y,
                               size_at_age,
                               return_data) {
    sp <- params@species_params
    sp <- set_species_param_default(sp, "age_mat", age_mat_vB(params))

    # size at age
    if (!is.null(size_at_age)) {
        if (!"species" %in% names(size_at_age)) {
            stop("The size at age data frame needs to have a 'species' column.")
        }
        if (!"age" %in% names(size_at_age)) {
            stop("The size at age data frame needs to have an 'age' column.")
        }
        if (!any(c("weight", "length") %in% names(size_at_age))) {
            stop("The size at age data frame needs to have either a 'length' or a 'weight' column.")
        }
        if (!"weight" %in% names(size_at_age)) {
            size_at_age <- left_join(size_at_age, select(sp, species, a, b),
                                     by = "species")
            size_at_age$weight <- size_at_age$a * size_at_age$length ^ size_at_age$b
        }
    }

    species <- valid_species_arg(params, species)
    # needed later to not confuse variable and column name
    selected_species <- species

    sp_sel <- sp$species %in% species
    ws <- getGrowthCurves(params, species, max_age, percentage)
    plot_dat <- reshape2::melt(ws)
    plot_dat$Species <- factor(plot_dat$Species, sp$species)
    plot_dat$Legend <- "model"

    # creating some VB
    if (all(c("a", "b", "k_vb", "w_inf") %in% names(sp))) {
        if ("t0" %in% names(sp)) {
            t0 <- sp$t0
        } else {
            t0 <- 0
        }
        VBdf <- data.frame("species" = sp$species,
                           "w_inf" = sp$w_inf,
                           "a" = sp$a,
                           "b" = sp$b,
                           "k_vb" = sp$k_vb,
                           "t0" = t0)
        VBdf$L_inf <- (VBdf$w_inf / VBdf$a) ^ (1 / VBdf$b)
        plot_dat2 <- plot_dat
        plot_dat2$value <-
            apply(plot_dat, 1,
                  function(x) {
                      sel <- VBdf$species == x[1]
                      length <- VBdf$L_inf[sel] *
                          (1 - exp(-VBdf$k_vb[sel] *
                                       (as.numeric(x[2]) - VBdf$t0[sel])))
                      VBdf$a[sel] * length ^ VBdf$b[sel]
                  })
        plot_dat2$Legend <- "von Bertalanffy"
        plot_dat <- rbind(plot_dat, plot_dat2)
    }
    y_label <- if (percentage)
        "Percent of maximum size"
    else "Size [g]"
    names(plot_dat)[3] <- y_label
    if (return_data) return(plot_dat)

    legend_levels <-
        intersect(c(dimnames(params@initial_n)$sp,
                    "Background", "Resource", "Total"),
                  plot_dat$Species)
    plot_dat$Species <- factor(plot_dat$Species, levels = legend_levels)
    plot_dat$LineSpec <- plot_dat$Species
    linesize <- make_linesize(legend_levels, highlight)
    p <- ggplot(filter(plot_dat, Legend == "model"),
                aes(group = Species)) +
        geom_line(aes(x = Age, y = .data[[y_label]],
                      colour = .data[["LineSpec"]], linetype = .data[["LineSpec"]],
                      linewidth = .data[["LineSpec"]]))
    p <- p + scale_x_continuous(name = "Age [Years]",
                                trans = if (log_x) "log10" else "identity") +
        scale_y_continuous(name = y_label,
                           trans = if (log_y) "log10" else "identity") +
        scale_colour_manual(name = "Species", values = params@linecolour[legend_levels]) +
        scale_linetype_manual(name = "Species", values = params@linetype[legend_levels]) +
        scale_discrete_manual("linewidth", values = linesize, guide = "none")

    # starting cases now
    if (!percentage)  {
        if (length(species) == 1) {
            idx <- which(sp$species == species)
            w_mat <- sp$w_mat[idx]
            age_mat <- sp$age_mat[idx]
            p <- p + geom_hline(yintercept = w_mat, linetype = "dashed",
                                colour = "grey") +
                geom_vline(xintercept = age_mat, linetype = "dashed",
                           colour = "grey") +
                annotate("text", 0, w_mat, vjust = -1, label = "Maturity")
            if ("von Bertalanffy" %in% plot_dat$Legend) {
                p <- p + geom_line(data = filter(plot_dat,
                                                 Legend == "von Bertalanffy"),
                                   aes(x = Age, y = .data[[y_label]]))
            }
            if (!is.null(size_at_age)) {
                size_at_age <- filter(size_at_age, species == selected_species)
                p <- p + geom_point(aes(x = age, y = weight),
                                    data = size_at_age,
                                    alpha = 0.2)
            }

        } else if (species_panel) { # need to add either no panel if no param
            # for VB or create a panel without VB
            p <- ggplot(plot_dat) +
                geom_line(aes(x = Age, y = .data[[y_label]], colour = Legend)) +
                scale_x_continuous(name = "Age [years]",
                                   trans = if (log_x) "log10" else "identity") +
                scale_y_continuous(name = "Size [g]",
                                   trans = if (log_y) "log10" else "identity") +
                geom_hline(aes(yintercept = w_mat),
                           data = tibble(Species = factor(legend_levels),
                                         w_mat = sp$w_mat[sp_sel]),
                           linetype = "dashed",
                           colour = "grey") +
                geom_vline(aes(xintercept = age_mat),
                           data = tibble(Species = factor(legend_levels),
                                         age_mat = sp$age_mat[sp_sel]),
                           linetype = "dashed",
                           colour = "grey") +
                geom_hline(aes(yintercept = w_max),
                           data = tibble(Species = factor(legend_levels),
                                         w_max = sp$w_max[sp_sel]),
                           linetype = "solid",
                           colour = "grey") +
                facet_wrap(~Species, scales = "free_y")

        }

    }
    return(make_mizer_plot(p, c("Species", "Age", y_label)))

}

#' @rdname plotGrowthCurves
#' @usage NULL
#' @export
plotlyGrowthCurves <- function(object, species = NULL,
                               max_age = 20,
                               percentage = FALSE,
                               species_panel = FALSE,
                               highlight = NULL,
                               log_x = FALSE, log_y = FALSE, log = NULL) {
    argg <- as.list(environment())
    plotHover(do.call("plotGrowthCurves", argg))
}
