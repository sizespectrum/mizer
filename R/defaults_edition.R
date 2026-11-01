#' Default editions
#'
#' Function to set and get which edition of default choices is being used.
#'
#' The mizer functions for creating new models make a lot of choices for default
#' values for parameters that are not provided by the user. Sometimes we find
#' better ways to choose the defaults and update mizer accordingly. When we do
#' this, we will increase the edition number.
#'
#' If you call `defaults_edition()` without an argument it returns the
#' currently active edition. If you call it with a number it sets the active
#' edition to that value.
#'
#' Users who want their existing code for creating models not to change
#' behaviour when run with future versions of mizer should explicitly set the
#' desired defaults edition at the top of their code.
#'
#' @section Each model keeps its own edition:
#' The active edition is used when a new model is created, and the model
#' records it. Whenever mizer later fills in a default for that model, for
#' example when it recalculates the model after you have changed a species
#' parameter, it uses the edition recorded in the model, not the one that is
#' active at the time. So changing the active edition affects only models
#' created afterwards, and an existing model keeps behaving as it did when it
#' was created. Call `defaults_edition(params)` to see the edition of a
#' model. Models created with a version of mizer that did not yet record the
#' edition are of edition 1.
#'
#' To give an existing model the defaults of a different edition, create it
#' afresh with that edition active.
#'
#' `defaults_edition(params) <- value` changes the edition recorded in a
#' model. It does not recalculate anything: the model keeps the values it has,
#' and the new edition is used only for defaults that mizer fills in from then
#' on. It is meant for a model that was created with a different edition
#' active, before mizer recorded the edition in the model, and has therefore
#' been marked as edition 1.
#'
#' @section The editions:
#' The most recent edition is edition 2. It will become the default in the
#' next release. The current default is edition 1. The following defaults
#' are changed in edition 2:
#'
#' * `catchability` = 0.3 instead of 1
#' * `initial effort` = 1 instead of 0
#' * `gamma` is set to ensure a feeding level of `f0` for larvae with the
#'   current value of `interaction_resource` instead of
#'   interaction_resource = 1`.
#' * `initial_n` is set using [get_steady_state_n()] instead of the rather
#'   arbitrary old choice.
#' * In [setReproduction()], `psi` is no longer forced to 1 above
#'   `w_repro_max`; its value is determined entirely by the maturity ogive
#'   and the reproductive proportion.
#'
#' @param edition NULL, a number, or a MizerParams or MizerSim object.
#' @return If `edition` is `NULL`, the currently active edition number. If
#'   `edition` is a MizerParams or MizerSim object, the edition recorded in
#'   that model. If `edition` is a number, the function sets the global
#'   `mizer_defaults_edition` option, emits a message, and returns the supplied
#'   value invisibly.
#' @export
#' @examples
#' # The edition used for new models
#' defaults_edition()
#' # The edition of an existing model
#' defaults_edition(NS_params)
defaults_edition <- function(edition = NULL) {
    current_edition <- 1
    # if called with a model, return the edition recorded in it
    if (is(edition, "MizerSim")) {
        edition <- edition@params
    }
    if (is(edition, "MizerParams")) {
        return(model_defaults_edition(edition))
    }
    # if called without argument or with NULL, return the current edition
    if (is.null(edition)) {
        edition <- getOption("mizer_defaults_edition")
        return(ifelse(is.null(edition), current_edition, edition))
    }
    # else check validity and set option
    assert_that(is.numeric(edition),
                edition >= 1)
    options(mizer_defaults_edition = edition)
    message("Mizer parameter defaults are now at edition ", edition, ".")
    return(invisible(edition))
}

#' @rdname defaults_edition
#' @param params A MizerParams object.
#' @param value The edition to record in the model.
#' @return `defaults_edition<-`: The MizerParams object with the new edition
#'   recorded. Nothing is recalculated.
#' @export
`defaults_edition<-` <- function(params, value) {
    assert_that(is(params, "MizerParams"),
                is.number(value),
                value >= 1)
    params@defaults_edition <- value
    params@time_modified <- lubridate::now()
    params
}

# The defaults edition recorded in a model. A model created before mizer
# recorded the edition, and not yet upgraded, has none and is of edition 1.
model_defaults_edition <- function(params) {
    # Unlike a missing list element, a missing S4 slot errors when read.
    if (!.hasSlot(params, "defaults_edition")) return(1)
    edition <- params@defaults_edition
    if (is.null(edition)) 1 else edition
}

# Evaluates `code` with `edition` as the active defaults edition. Mizer changes
# some models by building a new model and copying the old one into it. That new
# model has to be built with the defaults of the old one, not with those of
# whatever edition happens to be active.
with_defaults_edition <- function(edition, code) {
    old <- options(mizer_defaults_edition = edition)
    on.exit(options(old), add = TRUE)
    code
}
