# Machinery for the information that mizer gives the user while it is setting
# up or changing a model. See `with_info_level()` for a description.

#' Collect and report the information signals raised while setting parameters
#'
#' While mizer sets up or changes a model it raises conditions of class
#' `info_about_default` to tell the user about the choices it made on their
#' behalf. Each condition carries a `var` naming the quantity it is about and a
#' `level` giving its importance, where a low level means important and a high
#' level means chatter. This function evaluates `expr` with a calling handler
#' that collects those conditions and reports them together once `expr` has
#' finished, so that the user gets one report rather than a stream of messages,
#' with one line per `var`.
#'
#' Conditions that also have class `info_about_frozen` are reported as a
#' **warning** rather than as a message. These say that a change the user made
#' had no effect on the model because the rate array it feeds has been set by
#' hand, see [signal_frozen()]. That is not chatter about a default but a
#' failed instruction, and it has to be a warning so that it survives the
#' `suppressMessages()` in [species_params<-()], which is there to quieten the
#' routine recalculation chatter.
#'
#' The `info_level` argument controls how much is reported:
#' \itemize{
#' \item A number: conditions with `level` at most `info_level` are reported,
#'   the rest are dropped. So `info_level = 0` is silence.
#' \item `NA`: nothing is reported and nothing is dropped. The conditions are
#'   passed on to a handler further out, which is what an inner call uses when
#'   its caller has installed its own handler and will do the reporting.
#' }
#'
#' @param expr The expression to evaluate. It is evaluated in the calling
#'   environment, so assignments made in it have the same effect as they would
#'   have without this wrapper.
#' @param info_level The level of information to report, or `NA` to leave the
#'   reporting to a handler further out.
#'
#' @return The value of `expr`.
#' @concept helper
with_info_level <- function(expr, info_level = 3) {
    infos <- list()
    warns <- list()
    collect_info <- function(cnd) {
        # With `NA` we neither collect nor muffle, so that the condition
        # reaches the handler installed by our caller.
        if (is.na(info_level)) {
            return()
        }
        if (cnd$level <= info_level) {
            if (inherits(cnd, "info_about_frozen")) {
                warns[[cnd$var]] <<- cnd$message
            } else {
                infos[[cnd$var]] <<- cnd$message
            }
        }
        # Muffle even the conditions that we do not report, because we have
        # taken responsibility for them: `info_level = 0` means silence.
        cnd_muffle(cnd)
    }
    result <- withCallingHandlers(expr, info_about_default = collect_info)
    if (length(infos) > 0) {
        message(paste(infos, collapse = "\n"))
    }
    if (length(warns) > 0) {
        warning(paste(warns, collapse = "\n"), call. = FALSE)
    }
    result
}

#' Signal that a change the user made cannot take effect
#'
#' A rate array that has been set by hand is protected by a comment, see the
#' "Setting or changing rates" section in [setParams()]. Mizer then no longer
#' calculates it from the species parameters, so a change to one of the species
#' parameters that feeds it has no effect on the model. This function raises
#' the condition that tells the user so.
#'
#' The condition has class `info_about_frozen`, which is a subclass of
#' `info_about_default`, so it is collected by the same handler as the
#' information about default values, see [with_info_level()], but is reported
#' as a warning rather than as a message. A warning is needed because the
#' documented way of changing a species parameter, [species_params<-()], runs
#' `suppressMessages()` over the recalculation.
#'
#' Only signal this when the user has actually asked for something that is not
#' happening. The mere fact that a frozen array differs from what the formula
#' would give is not enough: mizer freezes arrays itself when it builds the
#' trait-based and community models, and those arrays differ from the formula
#' for the lifetime of the model. See [signal_frozen_changes()], which decides
#' this from the species parameters the user changed.
#'
#' @param var A string naming the quantity the report is about. Reports are
#'   collected by `var`, so each quantity is reported only once per call.
#' @param message The message to give the user.
#'
#' @return `NULL` invisibly. Called for its side effect of signalling.
#' @concept helper
signal_frozen <- function(var, message) {
    inform(message,
           class = c("info_about_frozen", "info_about_default"),
           var = var, level = 1)
}

#' Signal that a rate array was not recalculated because it is frozen
#'
#' Raised by the rate setters when they leave a frozen array alone although the
#' species parameters say that it should have a different value. This is a
#' plain `info_about_default` condition, so it is reported as a message and
#' [with_info_level()] silences it along with the other information when
#' `info_level = 0`. Where no handler is installed, for example when a rate
#' setter is called directly rather than via [setParams()], it surfaces as an
#' ordinary message. The stronger [signal_frozen()] warning is raised
#' elsewhere, by whoever knows that the user asked for a change, see
#' [signal_frozen_changes()].
#'
#' @param var A string naming the slot that was not recalculated.
#' @param quantity A string naming the quantity for the user, for example
#'   "metabolic rate".
#' @param reset_call A string with the call that recalculates the quantity, for
#'   example "setMetabolicRate(params, reset = TRUE)".
#' @param derived_from A string naming the parameters that the quantity would
#'   have been calculated from.
#'
#' @return `NULL` invisibly. Called for its side effect of signalling.
#' @concept helper
signal_not_recalculated <- function(var, quantity, reset_call,
                                    derived_from = "species parameters") {
    inform(paste0(
        "The ", quantity, " has been set manually and so is not recalculated ",
        "from the ", derived_from, ". Call `", reset_call, "` if you want the ",
        quantity, " to be calculated from the ", derived_from, " again."),
        class = "info_about_default", var = var, level = 1)
}

#' Which species parameters feed which frozen rate array
#'
#' A lookup table used by [signal_frozen_changes()] to decide whether a change
#' the user made to the species parameters can take effect. Each entry is named
#' after a slot of \linkS4class{MizerParams} that can be frozen and gives the
#' quantity as the user knows it, the call that unfreezes it, and the species
#' parameters that the setter and its default calculations read.
#'
#' The list of parameters does not have to be exhaustive, and deliberately is
#' not: it names the parameters that the setter reads directly together with
#' the main inputs of the default calculations for those parameters. A
#' parameter that is missing simply means that the user is not warned, and is
#' left with the message that the setter itself gives, see
#' [signal_not_recalculated()]. Listing a parameter that in fact has no
#' influence is the worse mistake, because it warns about a change that did
#' take effect.
#'
#' @return A named list of lists with entries `quantity`, `reset_call` and
#'   `params`.
#' @concept helper
frozen_rate_params <- function() {
    list(
        intake_max = list(
            quantity = "maximum intake rate",
            reset_call = "setMaxIntakeRate(params, reset = TRUE)",
            params = c("h", "n")),
        search_vol = list(
            quantity = "search volume",
            reset_call = "setSearchVolume(params, reset = TRUE)",
            params = c("gamma", "q", "n", "f0", "h", "interaction_resource")),
        metab = list(
            quantity = "metabolic rate",
            reset_call = "setMetabolicRate(params, reset = TRUE)",
            params = c("ks", "k", "p", "n", "fc", "alpha")),
        mu_b = list(
            quantity = "external mortality rate",
            reset_call = "setExtMort(params, reset = TRUE)",
            params = c("z0", "z_ext", "d", "n", "w_inf", "l_inf")),
        ext_encounter = list(
            quantity = "external encounter rate",
            reset_call = "setExtEncounter(params, reset = TRUE)",
            params = c("E_ext", "n")),
        ext_diffusion = list(
            quantity = "external diffusion rate",
            reset_call = "setExtDiffusion(params, reset = TRUE)",
            params = c("D_ext", "n")),
        pred_kernel = list(
            quantity = "predation kernel",
            reset_call = "setPredKernel(params, reset = TRUE)",
            params = c("pred_kernel_type", "beta", "sigma", "kernel_exp",
                       "ppmr_min", "ppmr_max")),
        maturity = list(
            quantity = "maturity ogive",
            reset_call = "setReproduction(params, reset = TRUE)",
            params = c("w_mat", "w_mat25", "l_mat", "l_mat25")),
        psi = list(
            quantity = "reproductive proportion",
            reset_call = "setReproduction(params, reset = TRUE)",
            params = c("w_mat", "w_mat25", "l_mat", "l_mat25", "m", "n",
                       "w_repro_max", "l_repro_max", "w_max", "l_max"))
    )
}

#' Signal the changes to species parameters that cannot take effect
#'
#' Goes through the rate arrays that can be frozen, see [frozen_rate_params()],
#' and raises a [signal_frozen()] condition for each frozen array that one of
#' the changed species parameters feeds. This is what turns "the model no
#' longer follows the species parameters" into a warning the user sees at the
#' moment they make the change, see [species_params<-()].
#'
#' @param params A \linkS4class{MizerParams} object, holding the rate arrays as
#'   they are, that is, before the change is applied.
#' @param changed A character vector with the names of the species parameters
#'   that the user changed.
#'
#' @return `NULL` invisibly. Called for its side effect of signalling.
#' @concept helper
signal_frozen_changes <- function(params, changed) {
    if (length(changed) == 0) {
        return(invisible(NULL))
    }
    for (var in names(frozen_rate_params())) {
        entry <- frozen_rate_params()[[var]]
        if (is.null(comment(slot(params, var)))) {
            next
        }
        affected <- intersect(changed, entry$params)
        if (length(affected) == 0) {
            next
        }
        signal_frozen(var, paste0(
            "Your change to the species ",
            ngettext(length(affected), "parameter ", "parameters "),
            paste0("`", affected, "`", collapse = ", "),
            " has not taken effect because the ", entry$quantity,
            " has been set manually and so is no longer calculated from the ",
            "species parameters. Call `", entry$reset_call, "` if you want ",
            "the ", entry$quantity, " to be calculated from the species ",
            "parameters again."))
    }
    invisible(NULL)
}

#' Which species parameters changed
#'
#' Compares the columns of two species parameter data frames and returns the
#' names of those that differ, using the same entry-by-entry comparison as
#' [record_given_species_params()]: a column that is not present in `old_sp` is
#' new and therefore counts as changed, and `NA` is compared as a value rather
#' than as an unknown.
#'
#' @param value A data frame with the new species parameters.
#' @param old_sp A data frame with the species parameters as they were before.
#'
#' @return A character vector with the names of the changed columns.
#' @concept helper
changed_species_params <- function(value, old_sp) {
    no_sp <- nrow(value)
    cols <- names(value)
    changed <- vapply(cols, function(col) {
        any(changed_entries(value[[col]], old_sp[[col]], no_sp))
    }, logical(1))
    cols[changed]
}
