#' Set line colours to be used in mizer plots
#' 
#' Colours for names that already had a colour set will be overwritten by
#' the colour you specify. Colours for names that did not yet have a colour
#' will be appended to the list of colours.
#' @param params A MizerParams object
#' @param colours A named list or named vector of line colours.
#' 
#' @return `setColours`: The MizerParams object with updated line colours
#' @export
#' @examples
#' params <- setColours(NS_params, list("Resource" = "red", "Total" = "#0000ff"))
#' species_params(params)["Cod", "linecolour"] <- "black"
#' plotSpectra(params, total = TRUE)
#' getColours(params)
setColours <- function(params, colours) {
    assert_that(is(params, "MizerParams"))
    colours <- validColours(colours)
    params@linecolour <- unlist(
        modifyList(as.list(params@linecolour), colours))
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

#' Set linetypes to be used in mizer plots
#' 
#' Linetypes for names that already had a linetype set will be overwritten by
#' the linetype you specify. Linetypes for names that did not yet have a 
#' linetype will be appended to the list of linetypes.
#' @param params A MizerParams object
#' @param linetypes A named list or named vector of linetypes.
#' 
#' @return `setLinetypes()`: The MizerParams object with updated linetypes
#' @export
#' @examples
#' params <- setLinetypes(NS_params, list("Total" = "dotted"))
#' species_params(params)["Cod", "linetype"] <- "dashed"
#' plotSpectra(params, total = TRUE)
#' getLinetypes(params)
setLinetypes <- function(params, linetypes) {
    assert_that(is(params, "MizerParams"))
    linetypes <- validLinetypes(linetypes)
    params@linetype <- unlist(
        modifyList(as.list(params@linetype), as.list(linetypes)))
    params
}

#' @rdname setLinetypes
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
        warning("The following are not valid lineypes and will be ",
                "ignored: ",
                paste(linetypes[!valid], collapse = ", "),
                ". The valid values are: ",
                paste(list_of_types, collapse = ", "))
    }
    as.list(linetypes[valid])
}