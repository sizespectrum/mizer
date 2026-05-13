#' Retired defaults editions
#' 
#' `r lifecycle::badge("deprecated")`
#'
#' The `defaults_edition()` mechanism has been retired. Mizer now always uses
#' the newer defaults introduced for defaults edition 2. This function remains
#' temporarily as a deprecated compatibility shim for old scripts.
#' 
#' @details
#' Calls to `defaults_edition()` no longer change any behaviour. Instead,
#' model setup code that needs the old defaults should spell them out
#' explicitly.
#' 
#' The defaults that are now always used are:
#' 
#' * `catchability` = 0.3 instead of 1
#' * initial fishing effort = 1 instead of 0
#' * `gamma` is set to ensure a feeding level of `f0` for larvae with the
#'   current value of `interaction_resource`
#' * `initial_n` is set using [get_steady_state_n()] instead of the rather
#'   arbitrary old choice.
#' * In [setReproduction()], `psi` is no longer forced to 1 above
#'   `w_repro_max`; its value is determined entirely by the maturity ogive
#'   and the reproductive proportion.
#' 
#' To reproduce old model setup behaviour:
#'
#' * add `catchability = 1` explicitly to `gear_params` or to the gear columns
#'   in `species_params`
#' * set old zero fishing effort explicitly with `initial_effort(params) <- 0`
#'   and provide full effort vectors or arrays to [project()] when missing
#'   values should mean zero
#' * if the old `gamma` behaviour matters, calculate or supply `gamma`
#'   explicitly before setting non-default `interaction_resource`
#' * if the old initial abundances matter, set `initialN(params) <- ...`
#'   explicitly. The old formula was
#'   `N = N0 * w_max^(2 * n - q - 2 + a) * w^(-n - a)`, with densities outside
#'   the species size range set to zero
#' * if the old reproduction allocation matters, set a custom
#'   `repro_prop`/`psi` explicitly rather than relying on forcing above
#'   `w_repro_max`
#'
#' @param edition Ignored. It is only kept so old calls do not fail.
#' @return The number 2, invisibly if an edition was supplied.
#' @export
#' @concept deprecated
defaults_edition <- function(edition = NULL) {
    lifecycle::deprecate_warn(
        "3.0.0",
        "defaults_edition()",
        details = paste(
            "Mizer now always uses the newer defaults.",
            "Old model setup code should set any old default values explicitly."
        )
    )
    if (is.null(edition)) {
        return(2)
    }
    invisible(2)
}
