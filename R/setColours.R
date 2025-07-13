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
#' Do not use this for setting the colours or linetypes of species, because
#' those are determined by setting the `linecolour` and `linetype` variables in
#' the species parameter data frame.
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
#' # Set colours and linetypes for species
#' species_params(params)["Cod", "linecolour"] <- "black"
#' species_params(params)["Cod", "linetype"] <- "dashed"
#' plotSpectra(params, total = TRUE)
#' getColours(params)
#' getLinetypes(params)
setColours <- function(params, colours) {
    assert_that(is(params, "MizerParams"))
    colours <- validColours(colours)
    if (identical(colours, as.list(params@linecolour))) {
        return(params)
    }
    params@linecolour <- unlist(
        modifyList(as.list(params@linecolour), colours))
    
    params@time_modified <- lubridate::now()
    params
}

#' @rdname setColours
#' @return `getColours()`: A named vector of colours
#' @export
getColours <- function(params) {
    params@linecolour
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
    assert_that(is(params, "MizerParams"))
    linetypes <- validLinetypes(linetypes)
    if (identical(linetypes, as.list(params@linetype))) {
        return(params)
    }
    params@linetype <- unlist(
        modifyList(as.list(params@linetype), as.list(linetypes)))
    
    params@time_modified <- lubridate::now()
    params
}

#' @rdname setColours
#' @return `getLinetypes()`: A named vector of linetypes
#' @export
getLinetypes <- function(params) {
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
