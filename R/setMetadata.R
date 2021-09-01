#' Set metadata for a model
#' 
#' @section Setting model metadata:
#' 
#' The metadata of a model consists of the following slots:
#' * `title` A string with the title for the model
#' * `description` A string with a description of the model
#' * `authors` A vector of strings with the author names
#' * `orcid` A named vector of strings where each name is an author name and
#'    each value is an orcid id.
#' * `version` The version string of the mizer version that last modified the
#'    model as returned by packageVersion("mizer")
#' * `extensions` A named vector of strings where each name is the name of
#'    and extension package needed to run the model and each value is a string 
#'    giving the information that the remotes package needs to install the 
#'    correct version of the extension package, see https://remotes.r-lib.org/
#' * `time_created` A POSIXct date-time object with the creation time
#' * `time_modified` A POSIXct date-time object with the last modified time
#' 
#' The `version`, `time_created` `time_modified` and `extensions` fields are
#' set automatically by the appropriate functions. The other fields are set
#' via `setMetadata()`.
#'
#' @param params The MizerParams object for the model
#' @param title A string with the title for the model
#' @param description A string with a description of the model
#' @param authors A vector of strings with the author names
#' @param orcid A named vector of strings where each name is an author name and
#'    each value is an orcid id.
#' @return The MizerParams object with updated metadata
#' @export
setMetadata <- function(params, title = NULL, description = NULL,
                        authors = NULL, orcid = NULL, ...) {
    params <- validParams(params)
    if (!is.null(title)) {
        assert_that(is.string(title))
        params@title <- title
    }
    if (!is.null(description)) {
        assert_that(is.string(description))
        params@description <- description
    }
    if (!is.null(authors)) {
        assert_that(is.character(authors))
        params@authors <- authors
    }
    if (!is.null(orcid)) {
        assert_that(is.character(orcid))
        if (is.null(names(orcid))) {
            stop("`orcid` should be a named vector.")
        }
        params@orcid <- orcid
    }
    
    params@time_modified <- lubridate::now()
    params
}

#' @rdname setMetadata
#' @export
#' @return `getMetadata()`: A list with the entries `title`, `description`,
#'   `authors`, `orcid`, `version`, `extensions`, `time_created` and
#'   `time_modified`.
getMetadata <- function(params) {
    list(title = params@title,
         description = params@description,
         authors = params@authors,
         orcid = params@orcid,
         version = params@version,
         extensions = params@extensions,
         time_created = params@time_created,
         time_modified = params@time_modified)
}
