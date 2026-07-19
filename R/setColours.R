#' Set line colours and line types to be used in mizer plots
#'
#' `r lifecycle::badge("experimental")`
#' Used for setting the colour and type of lines representing "Total",
#' "Resource", "Fishing", "Background", "External" and possibly other categories
#' in plots.
#'
#' Colours for names that already had a colour set for them will be overwritten
#' by the colour you specify. Colours for names that did not yet have a colour
#' will be appended to the list of colours.
#'
#' If a name coincides with the name of a species, the `linecolour` (for
#' `setColours()`) or `linetype` (for `setLinetypes()`) entry for that species
#' in `species_params` and `given_species_params` is updated as well, so that
#' the choice persists with the species. Alternatively you can set the
#' `linecolour` and `linetype` variables in the species parameter data frame
#' directly, see the example below.
#'
#' You can use the same colours in your own ggplot2 plots by adding
#' `scale_colour_manual(values = getColours(params))` to your plot. Similarly
#' you can use the linetypes with
#' `scale_linetype_manual(values = getLinetypes(params))`.
#'
#' @param params A MizerParams object
#' @param colours A named list or named vector of line colours.
#'
#' @return `setColours`: The MizerParams object with updated line colours
#' @export
#' @examples
#' params <- setColours(NS_params, list("Resource" = "red","Total" = "#0000ff"))
#' params <- setLinetypes(NS_params, list("Total" = "dotted"))
#' # Set colours and linetypes for species, either via setColours()/
#' # setLinetypes() or directly via the species parameter data frame
#' params <- setColours(params, list("Cod" = "black"))
#' species_params(params)["Cod", "linetype"] <- "dashed"
#' plotSpectra(params, total = TRUE)
#' getColours(params)
#' getLinetypes(params)
setColours <- function(params, colours) {
    UseMethod("setColours")
}
#' @export
setColours.MizerParams <- function(params, colours) {
    colours <- validColours(colours)
    if (identical(colours, as.list(params@linecolour))) {
        return(params)
    }
    params@linecolour <- unlist(
        modifyList(as.list(params@linecolour), colours))

    params <- sync_species_plot_column(params, colours, "linecolour")

    params@time_modified <- lubridate::now()
    params
}

#' @rdname setColours
#' @return `getColours()`: A named vector of colours
#' @export
getColours <- function(params) {
    UseMethod("getColours")
}
#' @export
getColours.MizerParams <- function(params) {
    params@linecolour
}

# If any of the given names are species names, update the corresponding
# entry in `species_params` and `given_species_params` so that the choice
# persists with the species rather than living only in the plotting slot.
sync_species_plot_column <- function(params, values, column) {
    species_values <- values[names(values) %in% params@species_params$species]
    if (length(species_values) == 0) {
        return(params)
    }
    idx <- match(names(species_values), params@species_params$species)
    new_vals <- unlist(species_values)

    old_vals <- params@species_params[[column]]
    if (is.null(old_vals)) {
        old_vals <- rep(NA_character_, nrow(params@species_params))
    }
    old_vals[idx] <- new_vals
    params@species_params[[column]] <- old_vals

    old_given <- params@given_species_params[[column]]
    if (is.null(old_given)) {
        old_given <- rep(NA_character_, nrow(params@given_species_params))
    }
    old_given[idx] <- new_vals
    params@given_species_params[[column]] <- old_given

    params
}

validColours <- function(colours) {
    valid <- sapply(colours, function(X) {
        tryCatch(is.matrix(col2rgb(X)),
                 error = function(e) FALSE)
    })
    if (!all(valid)) {
        warning("The following are not valid colour values and will be ",
                "ignored: ",
                paste(colours[!valid], collapse = ", "))
    }
    as.list(colours[valid & !is.na(colours)])
}


#' @rdname setColours
#' @param linetypes A named list or named vector of linetypes.
#'
#' @return `setLinetypes()`: The MizerParams object with updated linetypes
#' @export
setLinetypes <- function(params, linetypes) {
    UseMethod("setLinetypes")
}
#' @export
setLinetypes.MizerParams <- function(params, linetypes) {
    linetypes <- validLinetypes(linetypes)
    if (identical(linetypes, as.list(params@linetype))) {
        return(params)
    }
    params@linetype <- unlist(
        modifyList(as.list(params@linetype), as.list(linetypes)))

    params <- sync_species_plot_column(params, linetypes, "linetype")

    params@time_modified <- lubridate::now()
    params
}

#' @rdname setColours
#' @return `getLinetypes()`: A named vector of linetypes
#' @export
getLinetypes <- function(params) {
    UseMethod("getLinetypes")
}
#' @export
getLinetypes.MizerParams <- function(params) {
    params@linetype
}



validLinetypes <- function(linetypes) {
    linetypes <- linetypes[!is.na(linetypes)]
    list_of_types <- list(0, 1, 2, 3, 4, 5, 6, "blank", "solid", "dashed",
                          "dotted", "dotdash", "longdash", "twodash")
    valid <- linetypes %in% list_of_types

    if (!all(valid)) {
        warning("The following are not valid linetypes and will be ",
                "ignored: ",
                paste(linetypes[!valid], collapse = ", "),
                ". The valid values are: ",
                paste(list_of_types, collapse = ", "))
    }
    as.list(linetypes[valid])
}
