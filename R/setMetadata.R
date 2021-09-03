#' Set metadata for a model
#' 
#' Setting metadata is particularly important for sharing your model with
#' others. All metadata fields are optional and you can also add other fields
#' of your own choosing. If you set a value
#' for a field that already existed, the old value will be overwritten.
#' 
#' In addition to the metadata fields you can set by hand, there are four fields
#' that are set automatically by mizer:
#' 
#' * `mizer_version` The version string of the mizer version under which the
#' model was last saved. Can be compared to the current version which is
#' obtained with `packageVersion("mizer")`. The purpose of this field is that
#' if the model is not working as expected in the current version of mizer,
#' you can go back to the older version under which presumably it was working.
#' * `extensions` A named vector of strings where each name is the name of and
#' extension package needed to run the model and each value is a string giving
#' the information that the remotes package needs to install the correct version
#' of the extension package, see https://remotes.r-lib.org/. This field is
#' set by the extension packages.
#' * `time_created` A POSIXct date-time object with the creation time.
#' * `time_modified` A POSIXct date-time object with the last modified time.
#'
#' @param params The MizerParams object for the model
#' @param title A string with the title for the model
#' @param description A string with a description of the model. This could for
#'   example contain information about any publications using the model.
#' @param authors An author entry or a list of author entries, where each author
#'   entry could either be just a name or could itself be a list with fields
#'   like `name`, `orcid`, possibly `email`.
#' @param url A URL where more information about the model can be found. This
#'    could be a blog post on the mizer blog, for example.
#' @param doi The digital object identifier for your model. To create a doi you
#'   can use online services like https://zenodo.org/ or https://figshare.com.
#' @param ... Additional metadata fields that you would like to add
#' 
#' @return `setMetadata()`: The MizerParams object with updated metadata
#' @export
setMetadata <- function(params, title, description,
                        authors, url, doi, ...) {
    params <- validParams(params)
    extra <- list(...)
    special <- c("mizer_version", "extensions", "time_modified", "time_created")
    if (any(special %in% names(extra))) {
        message("The fields ", special, " are set automatically by mizer. ",
                "The values you supplied will be ignored.")
        extra <- extra[!(names(extra) %in% special)]
    }
    
    metadata <- params@metadata
    if (!missing(title)) {
        assert_that(is.string(title))
        metadata$title <- title
    }
    if (!missing(description)) {
        assert_that(is.string(description))
        metadata$description <- description
    }
    if (!missing(authors)) {
        metadata$authors <- authors
    }
    if (!missing(url)) {
        assert_that(is.character(url))
        pmetadata$url <- url
    }
    if (!missing(doi)) {
        assert_that(is.character(doi))
        metadata$doi <- doi
    }
    params@metadata <- modifyList(metadata, list(...))
    
    params@time_modified <- lubridate::now()
    params
}

#' @rdname setMetadata
#' @export
#' @return `getMetadata()`: A list with all metadata entries that have been set,
#'   including at least 
#'   `mizer_version`, `extensions`, `time_created` and `time_modified`.
getMetadata <- function(params) {
    list <- params@metadata
    list$mizer_version = params@mizer_version
    list$extensions = params@extensions
    list$time_created = params@time_created
    list$time_modified = params@time_modified
    list
}
