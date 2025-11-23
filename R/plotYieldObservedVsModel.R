#' Plotting observed vs. model yields
#'
#' `r lifecycle::badge("experimental")`
#' If yield observations are available for at least some species via the
#' `yield_observed` column in the species parameter data frame, this function
#' plots the yield of each species in the model against the observed
#' yields. When called with a MizerSim object, the plot will use the model
#' yields predicted for the final time step in the simulation.
#'
#' Before you can use this function you will need to have added a
#' `yield_observed` column to your model which gives the observed yield in
#' grams per year.  For species for which you have no observed yield, you should set
#' the value in the `yield_observed` column to 0 or NA.
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
#' @param ... Additional arguments passed to the generic function.
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
#' species_params(params)$yield_observed <-
#'     c(0.8, 61, 12, 35, 1.6, NA, 10, 7.6, 135, 60, 30, NA)
#' params <- calibrateYield(params)
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
    
    biomass <- sweep(params@initial_n, 2, params@w * params@dw, "*")
    yield_model <- rowSums(biomass * getFMort(params))
    
    # Select appropriate species
    species <- valid_species_arg(object, species)
    no_yield <- yield_model[species] == 0
    if (any(no_yield)) {
        message("The following species are not being fished in your model ",
                "and will not be included in the plot: ",
                paste0(species[no_yield], collapse = ", "), ".")
        species <- species[!no_yield]
    }
    if (length(species) == 0) stop("No species selected, please fix.")
    
    # find rows corresponding to species selected
    row_select <- match(species, sp_params$species)
    if (!"yield_observed" %in% names(sp_params)) {
        stop("You have not provided values for the column 'yield_observed' ",
             "in the mizerParams/mizerSim object.")
    } else if (!is.numeric(sp_params$yield_observed)) {
        stop("The column 'yield_observed' in the mizerParams/mizerSim object",
             " is not numeric, please fix.")
    } else { # accept
        yield_observed <- sp_params$yield_observed
    }
    
    # Build dataframe
    dummy <- data.frame(species = species, 
                       model = yield_model[row_select], 
                       observed = yield_observed[row_select]) %>%
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
    gg
}
#' @export
plotYieldObservedVsModel.MizerSim <- function(object, species = NULL, ratio = FALSE,
                                     log_scale = TRUE, return_data = FALSE, 
                                     labels = TRUE, show_unobserved = FALSE, ...) {
    params <- setInitialValues(object@params, object)
    plotYieldObservedVsModel(params, species = species, ratio = ratio,
                             log_scale = log_scale, return_data = return_data,
                             labels = labels, show_unobserved = show_unobserved)
}


#' @rdname plotYieldObservedVsModel
#' @export
plotlyYieldObservedVsModel <- function(object, species = NULL, ratio = FALSE,
                                         log_scale = TRUE, return_data = FALSE, 
                                         show_unobserved = FALSE, ...) {
    argg <- as.list(environment())
    argg$labels <- FALSE
    ggplotly(do.call("plotYieldObservedVsModel", argg), ...)
}