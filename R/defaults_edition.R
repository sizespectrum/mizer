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
#' currently active edition. Otherwise it sets the active edition to the 
#' given value.
#' 
#' Users who want their existing code for creating models not to change 
#' behaviour when run with future versions of mizer should explicitly set the
#' desired defaults edition at the top of their code.
#' 
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
#' @param edition NULL or a numerical value.
#' @return If `edition` is `NULL`, the currently active edition number. If
#'   `edition` is supplied, the function sets the global
#'   `mizer_defaults_edition` option, emits a message, and returns the supplied
#'   value invisibly.
#' @export
defaults_edition <- function(edition = NULL) {
    current_edition <- 1
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
