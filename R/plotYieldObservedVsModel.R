#' Plotting observed vs. model yields
#'
#' `r lifecycle::badge("experimental")`
#' If yield observations are available for at least some species via the
#' `yield_observed` column, this function plots the yield of each species in
#' the model against the observed yields. When called with a MizerSim object,
#' the plot will use the model yields predicted for the final time step in the
#' simulation.
#'
#' Before you can use this function you will need to have added a
#' `yield_observed` column to your model which gives the observed yield in
#' grams per year. Its home is the gear parameter data frame, see
#' [gear_params()], where you give the yield for each gear-species pair and
#' this function adds them up over the gears. For backwards compatibility a
#' `yield_observed` column in the species parameter data frame is also
#' accepted, see [get_yield_observed()]. For species for which you have no
#' observed yield, you should set the value in the `yield_observed` column to
#' 0 or NA.
#'
#' The total relative error is shown in the caption of the plot, calculated by
#'   \deqn{TRE = \sum_i|1-\rm{ratio_i}|}{TRE = sum_i |1-ratio_i|} where
#'   \eqn{\rm{ratio_i}}{ratio_i} is the ratio of model yield / observed
#'   yield for species i.
#'
#' @param object An object of class \linkS4class{MizerParams} or
#'   \linkS4class{MizerSim}.
#' @param species The species to be included. Optional. By default all observed
#'   yields will be included. A vector of species names, or a numeric vector
#'   with the species indices, or a logical vector indicating for each species
#'   whether it is to be included (TRUE) or not.
#' @param ratio Whether to plot model yield vs. observed yield (FALSE) or
#'   the ratio of model : observed yield (TRUE). Default is FALSE.
#' @param log_scale Whether to plot on the log10 scale (TRUE) or not (FALSE).
#'   For the non-ratio plot this applies for both axes, for the ratio plot only
#'   the x-axis is on the log10 scale. Default is TRUE.
#' @param labels Whether to show text labels for each species (TRUE) or not
#'   (FALSE). Default is TRUE.
#' @param show_unobserved Whether to include also species for which no
#'   yield observation is available. If TRUE, these species will be 
#'   shown as if their observed yield was equal to the model yield.
#' @param return_data Whether to return the data frame for the plot (TRUE) or
#'   not (FALSE). Default is FALSE.
#' @param ... For [plotlyYieldObservedVsModel()], additional arguments passed
#'   to [plotHover()]. Otherwise unused.
#' 
#' @return A ggplot2 object with the plot of model yield by species compared
#'   to observed yield. If `return_data = TRUE`, the data frame used to
#'   create the plot is returned instead of the plot.
#' @importFrom stats cor.test
#' @importFrom utils data
#' @export
#' @examples
#' # create an example
#' params <- NS_params
#' # In this model each species is caught by a single gear, so there is one
#' # row in the gear parameters for each species, in the same order.
#' # Species without an observation get NA.
#' gear_params(params)$yield_observed <-
#'     c(NA, NA, NA, 3e11, 4e9, 4e10, 5e10, NA, 2e11, 6e10, 3e11, NA)
#'
#' # Plot with default options
#' plotYieldObservedVsModel(params)
#' 
#' # Plot including also species without observations
#' plotYieldObservedVsModel(params, show_unobserved = TRUE)
#'
#' # Show the ratio instead
#' plotYieldObservedVsModel(params, ratio = TRUE)
plotYieldObservedVsModel <- function(object, species = NULL, ratio = FALSE,
                                     log_scale = TRUE, return_data = FALSE, 
                                     labels = TRUE, show_unobserved = FALSE, ...) {
    UseMethod("plotYieldObservedVsModel")
}
#' @export
plotYieldObservedVsModel.MizerParams <- function(object, species = NULL, ratio = FALSE,
                                     log_scale = TRUE, return_data = FALSE, 
                                     labels = TRUE, show_unobserved = FALSE, ...) {
    
    params <- object
    sp_params <- params@species_params
    
    yield_model <- getYield(params)

    # Select appropriate species
    species <- valid_species_arg(object, species)
    no_yield <- yield_model[species] == 0
    if (any(no_yield)) {
        signal_info("yield_model", paste0(
            "The following species are not being fished in your model and ",
            "will not be included in the plot: ",
            paste0(species[no_yield], collapse = ", "), "."),
            level = 1, unhandled = "show")
        species <- species[!no_yield]
    }
    if (length(species) == 0) stop("No species selected, please fix.")
    
    # find rows corresponding to species selected
    row_select <- match(species, sp_params$species)
    yield_observed <- get_yield_observed(params)
    if (is.null(yield_observed)) {
        stop("You have not provided values for the column 'yield_observed' ",
             "in either the gear parameters or the species parameters of the ",
             "mizerParams/mizerSim object.")
    }

    # Build dataframe
    dummy <- data.frame(species = species,
                       model = yield_model[row_select],
                       observed = unname(yield_observed[row_select])) %>%
        mutate(species = factor(species, levels = species),
               is_observed = !is.na(observed) & observed > 0,
               observed = case_when(is_observed ~ observed, !is_observed ~ model),
               ratio = model/observed)
    
    
    # Check that at least one observed yield exists
    if (sum(dummy$is_observed) == 0) {
        warning("There are no observed yields to compare to model, ", 
                  "only plotting model yields.")
    }
    
    if (!show_unobserved) {
        dummy <- filter(dummy, is_observed)
    }
    
    if (return_data == TRUE) return(dummy)
    
    # Calculate total sum of differences (abs(1-ratio)) rounded to 3 digits
    tre <- round(sum(abs(1 - dummy$ratio)), digits = 3)
    
    caption <- paste0("Total relative error = ", tre)
    if (any(!dummy$is_observed)) {
        caption <- paste(caption,
                         "\n Open circles represent species without yield observation.")
    }
    
    if (ratio == FALSE) {
        gg <- ggplot(data = dummy, aes(x = observed, y = model,
                                       colour = species, shape = is_observed)) + 
            geom_abline(aes(intercept = 0, slope = 1), colour = 'purple',
                        linetype = "dashed", linewidth = 1.3) + # y = x line
            geom_point(size = 3) +
            labs(y = 'model yield [g/year]') +
            coord_cartesian(ylim = range(dummy$model, dummy$observed))
    } else {
        gg <- ggplot(data = dummy, aes(x = observed, y = ratio,
                                       colour = species, shape = is_observed)) + 
            geom_hline(aes(yintercept = 1), linetype = "dashed",
                       colour = 'purple', linewidth = 1.3) +
            geom_point(size = 3) +
            labs(y = 'model yield / observed yield') +
            coord_cartesian(ylim = range(dummy$ratio))
    }
    
    gg <- gg + labs(x = 'observed yield [g/year]', caption = caption) +
        scale_colour_manual(values = getColours(params)[as.character(dummy$species)]) +
        scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1)) +
        guides(shape = "none")
    
    
    if (log_scale == TRUE & ratio == FALSE) {
        gg = gg + scale_x_log10() + scale_y_log10()
    }
    if (log_scale == TRUE & ratio == TRUE) {
        gg = gg + scale_x_log10()
    }
    
    if (labels == TRUE)  {
        gg = gg + ggrepel::geom_label_repel(
            aes(label = species),
            box.padding = 0.35,
            point.padding = 0.5,
            segment.color = 'grey50',
            show.legend = FALSE,
            max.overlaps = Inf,
            seed = 42)
    }
    make_mizer_plot(gg, "all")
}
#' @export
plotYieldObservedVsModel.MizerSim <- function(object, species = NULL, ratio = FALSE,
                                     log_scale = TRUE, return_data = FALSE, 
                                     labels = TRUE, show_unobserved = FALSE, ...) {
    params <- finalParams(object)
    plotYieldObservedVsModel(params, species = species, ratio = ratio,
                             log_scale = log_scale, return_data = return_data,
                             labels = labels, show_unobserved = show_unobserved)
}


#' @rdname plotYieldObservedVsModel
#' @usage NULL
#' @export
plotlyYieldObservedVsModel <- function(object, species = NULL, ratio = FALSE,
                                         log_scale = TRUE,
                                         show_unobserved = FALSE, ...) {
    argg <- as.list(environment())
    argg$labels <- FALSE
    plotHover(do.call("plotYieldObservedVsModel", argg), ...)
}


#' Observed yield of each species
#'
#' The observed yield lives in the `yield_observed` column of the gear
#' parameter data frame, see [gear_params()], where it is given for each
#' gear-species pair. This function adds the observations up over the gears to
#' give the total observed yield of each species.
#'
#' Older models, and the examples in older versions of mizer, put
#' `yield_observed` into the species parameter data frame instead. That is
#' still accepted: a species that has no observation among the gear parameters
#' takes its value from the species parameters. Where both tables give a value
#' for a species, the gear parameters win.
#'
#' @param params A MizerParams object
#' @return A numeric vector with one entry for each species, named by species,
#'   holding the observed yield in grams per year, or `NA` for species without
#'   an observation. `NULL` if neither the gear parameters nor the species
#'   parameters have a `yield_observed` column.
#' @concept helper
get_yield_observed <- function(params) {
    species <- as.character(params@species_params$species)
    gp <- params@gear_params
    sp <- params@species_params
    in_gear <- "yield_observed" %in% names(gp)
    in_species <- "yield_observed" %in% names(sp)
    if (!in_gear && !in_species) return(NULL)

    yield <- rep(NA_real_, length(species))
    names(yield) <- species

    if (in_gear) {
        observed <- gp$yield_observed
        if (!is.numeric(observed)) {
            stop("The column 'yield_observed' in the gear parameter data ",
                 "frame is not numeric, please fix.")
        }
        given <- !is.na(observed)
        if (any(given)) {
            # Sum over all gears catching a species. Species without any
            # observation get NA from `tapply()`.
            totals <- tapply(observed[given],
                             factor(as.character(gp$species)[given],
                                    levels = species),
                             sum)
            yield[] <- as.numeric(totals)
        }
    }

    if (in_species) {
        observed <- sp$yield_observed
        if (!is.numeric(observed)) {
            stop("The column 'yield_observed' in the species parameter data ",
                 "frame is not numeric, please fix.")
        }
        no_obs <- is.na(yield)
        yield[no_obs] <- observed[no_obs]
    }

    yield
}
