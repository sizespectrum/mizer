# Machinery for the information that mizer gives the user while it is setting
# up or changing a model. See `with_info_level()` for a description.

# Package-local state recording whether a `with_info_level()` handler is
# already collecting. It is what makes the handlers nest: the outermost one
# does the reporting and the inner ones step aside, so that a function can wrap
# its body without having to know whether its caller has done the same. It is
# always restored on exit, including when the expression throws.
info_reporting <- new.env(parent = emptyenv())
info_reporting$active <- FALSE

#' The default level of information that mizer gives
#'
#' Returns the `mizer_info_level` option if it is set and `fallback`
#' otherwise. This is the default of the `info_level` argument of the functions
#' that report information, so that
#' `options(mizer_info_level = 0)` quietens mizer as a whole, including the
#' functions that have no `info_level` argument of their own, such as
#' [species_params<-()] and the rate setters.
#'
#' Extension packages should use this as the default of their own `info_level`
#' argument, so that a constructor or setter of theirs follows the option like
#' mizer's own do:
#'
#' ```
#' newFooParams <- function(species_params, ...,
#'                          info_level = default_info_level()) {
#'     newMultispeciesParams(species_params, info_level = info_level, ...)
#' }
#' ```
#'
#' Take the argument explicitly like that rather than hard-coding a value in the
#' call, which would make a user's own `info_level` collide with it.
#'
#' @param fallback The level to use when the option is not set. Defaults to 3,
#'   which reports everything.
#'
#' @return A single number, or `NA` to leave the reporting to a handler further
#'   out.
#' @export
#' @concept helper
#' @examples
#' default_info_level()
#'
#' # Setting the option changes what every reporting function defaults to.
#' old <- options(mizer_info_level = 1)
#' default_info_level()
#' options(old)
default_info_level <- function(fallback = 3) {
    getOption("mizer_info_level", default = fallback)
}

#' Collect and report the information signals raised while setting parameters
#'
#' While mizer sets up or changes a model it raises conditions of class
#' `info_about_default` to tell the user about the choices it made on their
#' behalf and about the instructions it could not carry out. This function
#' evaluates `expr` with a calling handler that collects those conditions and
#' reports them together once `expr` has finished, so that the user gets one
#' report rather than a stream of messages.
#'
#' Each condition carries three fields, see [signal_info()]:
#' \itemize{
#' \item `var` names the quantity the report is about.
#' \item `level` says how important it is, a low level meaning important and a
#'   high level meaning chatter. Only conditions with `level` at most
#'   `info_level` are reported.
#' \item `severity` says how to report it: `"info"` conditions become a single
#'   `message()` and `"warning"` conditions a single `warning()`.
#' }
#'
#' The severity matters because [species_params<-()] runs `suppressMessages()`
#' over its recalculation to quieten the routine chatter. A report that the
#' user needs to see even there — that an instruction of theirs had no effect —
#' must therefore be a warning, see [signal_frozen()].
#'
#' Identical reports are collapsed, so a quantity that is reported on twice in
#' the same call takes up one line, but two different things said about the
#' same quantity are both kept.
#'
#' # Nesting
#'
#' Handlers nest by themselves: while one is collecting, any handler installed
#' further in steps aside and lets the outer one do the reporting. A function
#' can therefore wrap its body in `with_info_level()` without knowing whether
#' its caller has already done so, which is what allows every entry point to
#' install a handler. `info_level = NA` asks for the same thing explicitly,
#' for the rare case where a function wants to leave the reporting to a caller
#' that has not installed a handler yet.
#'
#' The reporting happens on exit, so a function can wrap its whole body even
#' though it returns from the middle of it.
#'
#' Silence is the exception to "the outermost handler decides":
#' `info_level = 0` drops the reports raised inside it even when a handler
#' further out is collecting, so that a function can build something quietly
#' as part of a larger job that does report.
#'
#' @param expr The expression to evaluate. It is evaluated in the calling
#'   environment, so assignments made in it have the same effect as they would
#'   have without this wrapper.
#' @param info_level The level of information to report, or `NA` to leave the
#'   reporting to a handler further out. Defaults to
#'   [default_info_level()], which consults the `mizer_info_level` option.
#' @param except A character vector of `var`s not to report on, for a caller
#'   that is going to say the same thing itself. Everything else raised inside
#'   `expr` is reported as usual.
#'
#' @return The value of `expr`.
#' @export
#' @concept helper
#' @examples
#' # Wrap the body of a function that reports, and everything raised inside it
#' # is collected and given together once the call has finished.
#' myConstructor <- function(x, info_level = default_info_level()) {
#'     with_info_level(info_level = info_level, {
#'         signal_info("h", "No `h` provided, using a default.", level = 1)
#'         signal_info("gamma", "Calculating `gamma` from `f0`.")
#'         x
#'     })
#' }
#' myConstructor(1)
#'
#' # `info_level = 1` keeps only the report that was marked important.
#' myConstructor(1, info_level = 1)
#'
#' # `info_level = 0` is silence.
#' myConstructor(1, info_level = 0)
with_info_level <- function(expr, info_level = default_info_level(),
                            except = character()) {
    # `info_level = 0` silences its expression whatever a handler further out
    # would have done with the reports, so that a function can build something
    # quietly inside one that is reporting.
    silent <- isTRUE(info_level == 0)
    # Otherwise, if a handler further out is already collecting, or we were
    # told to leave the reporting to one, evaluate without collecting or
    # muffling anything. A caller that wants to drop a particular report has
    # to install a handler for it even so.
    if (!silent && length(except) == 0 &&
        (info_reporting$active ||
         !is.numeric(info_level) || length(info_level) != 1 ||
         is.na(info_level))) {
        return(expr)
    }
    reports <- list()
    collect_info <- function(cnd) {
        # The fields are defaulted rather than required, so that a condition
        # raised by an extension package that does not use `signal_info()` is
        # still reported rather than throwing.
        level <- if (is.null(cnd$level)) 3 else cnd$level
        severity <- if (is.null(cnd$severity)) "info" else cnd$severity
        if (level <= info_level && !isTRUE(cnd$var %in% except)) {
            key <- paste(severity, cnd$var, cnd$message, sep = "\r")
            reports[[key]] <<- list(severity = severity,
                                    message = cnd$message)
        }
        # Muffle even the conditions that we do not report, because we have
        # taken responsibility for them: `info_level = 0` means silence.
        cnd_muffle(cnd)
    }
    # The reporting is done on exit so that it also happens when `expr` leaves
    # early. A function can then wrap its whole body in `with_info_level()`
    # even though it returns from the middle of it, which is what makes this
    # usable as a wrapper around an existing function.
    was_active <- info_reporting$active
    info_reporting$active <- TRUE
    on.exit({
        info_reporting$active <- was_active
        severities <- vapply(reports, `[[`, character(1), "severity")
        messages <- vapply(reports, `[[`, character(1), "message")
        if (any(severities == "info")) {
            message(paste(messages[severities == "info"], collapse = "\n"))
        }
        if (any(severities == "warning")) {
            warning(paste(messages[severities == "warning"], collapse = "\n"),
                    call. = FALSE)
        }
    }, add = TRUE)
    withCallingHandlers(expr, info_about_default = collect_info)
}

#' Signal information about a choice mizer made
#'
#' Raises the condition that [with_info_level()] collects. This is the way for
#' mizer, and for anything extending it, to tell the user about a default it
#' filled in, an input it adjusted or an instruction it could not carry out,
#' without deciding on its own how loudly to say it: the handler installed by
#' whichever function the user actually called does that.
#'
#' Progress reports are the one thing that does not belong here: they have to
#' appear while the work is going on, and these are collected and given at the
#' end.
#'
#' @param var A string naming the quantity the report is about.
#' @param message The message to give the user.
#' @param level How important the report is. Level 1 is important enough to
#'   survive `info_level = 1`, level 3 is chatter that only the default
#'   `info_level = 3` shows.
#' @param severity `"info"` to report as a message, `"warning"` to report as a
#'   warning. Use `"warning"` when the user asked for something that is not
#'   happening, because a message can be, and on the
#'   [species_params<-()] path is, suppressed.
#' @param unhandled What to do when no handler is collecting, for example
#'   because a rate setter was called directly rather than through
#'   [setParams()]. `"drop"` says nothing, which suits chatter that only makes
#'   sense as part of a report about a whole model. `"show"` reports it there
#'   and then, at the same severity: a message for `"info"` and a warning for
#'   `"warning"`.
#' @param class Further classes to give the condition, for code that wants to
#'   catch a particular kind of report.
#'
#' @return `NULL` invisibly. Called for its side effect of signalling.
#' @export
#' @concept helper
#' @examples
#' # With nothing collecting, a `"drop"` report says nothing at all ...
#' signal_info("h", "Using a default for `h`.")
#'
#' # ... whereas `unhandled = "show"` reports it there and then.
#' signal_info("h", "Using a default for `h`.", unhandled = "show")
#'
#' # Normally it is raised inside a call whose body is wrapped in
#' # `with_info_level()`, which is what decides whether to show it.
#' with_info_level(signal_info("h", "Using a default for `h`."))
signal_info <- function(var, message, level = 3,
                        severity = c("info", "warning"),
                        unhandled = c("drop", "show"),
                        class = character()) {
    severity <- match.arg(severity)
    unhandled <- match.arg(unhandled)
    emit <- if (unhandled == "drop") {
        signal
    } else if (severity == "warning") {
        rlang::warn
    } else {
        inform
    }
    emit(message, class = c(class, "info_about_default"),
         var = var, level = level, severity = severity)
}

#' Signal that a change the user made cannot take effect
#'
#' A rate array that has been set by hand is protected by a comment, see the
#' "Setting or changing rates" section in [setParams()]. Mizer then no longer
#' calculates it from the species parameters, so a change to one of the species
#' parameters that feeds it has no effect on the model. This function raises
#' the condition that tells the user so.
#'
#' The condition is raised at severity `"warning"`, see [signal_info()], so
#' that it survives the `suppressMessages()` that [species_params<-()] runs
#' over its recalculation. It also carries the class `info_about_frozen` for
#' code that wants to catch this kind of report in particular.
#'
#' Only signal this when the user has actually asked for something that is not
#' happening. The mere fact that a frozen array differs from what the formula
#' would give is not enough: mizer freezes arrays itself when it builds the
#' trait-based and community models, and those arrays differ from the formula
#' for the lifetime of the model. See [signal_frozen_changes()], which decides
#' this from the species parameters the user changed.
#'
#' @param var A string naming the quantity the report is about.
#' @param message The message to give the user.
#'
#' @return `NULL` invisibly. Called for its side effect of signalling.
#' @concept helper
signal_frozen <- function(var, message) {
    signal_info(var, message, level = 1, severity = "warning",
                unhandled = "show", class = "info_about_frozen")
}

#' Signal that a rate array was not recalculated because it is frozen
#'
#' Raised by the rate setters when they leave a frozen array alone although the
#' species parameters say that it should have a different value. It is reported
#' as a message, and `info_level = 0` silences it along with the other
#' information. Where no handler is collecting, for example when a rate setter
#' is called directly rather than via [setParams()], it is shown anyway,
#' because it may then be all the user hears. The stronger [signal_frozen()]
#' warning is raised elsewhere, by whoever knows that the user asked for a
#' change, see [signal_frozen_changes()].
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
#' @export
#' @concept helper
#' @examples
#' with_info_level(
#'     signal_not_recalculated("metab", "metabolic rate",
#'                             "setMetabolicRate(params, reset = TRUE)")
#' )
signal_not_recalculated <- function(var, quantity, reset_call,
                                    derived_from = "species parameters") {
    signal_info(var, paste0(
        "The ", quantity, " has been set manually and so is not recalculated ",
        "from the ", derived_from, ". Call `", reset_call, "` if you want the ",
        quantity, " to be calculated from the ", derived_from, " again."),
        level = 1, unhandled = "show")
}

#' Which parameters feed which frozen array
#'
#' A lookup table used by [signal_frozen_changes()] to decide whether a change
#' the user made can take effect. Each entry is named after a slot of
#' \linkS4class{MizerParams} that can be frozen and gives the quantity as the
#' user knows it, the call that unfreezes it, the parameters that the setter
#' and its default calculations read, and what kind of parameters those are.
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
#' @return A named list of lists with entries `quantity`, `reset_call`,
#'   `params` and `derived_from`.
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
                       "w_repro_max", "l_repro_max", "w_max", "l_max")),
        cc_pp = list(
            quantity = "resource capacity",
            reset_call = "setResource(params, reset = TRUE)",
            params = c("kappa", "lambda", "w_pp_cutoff"),
            derived_from = "resource parameters"),
        rr_pp = list(
            quantity = "resource rate",
            reset_call = "setResource(params, reset = TRUE)",
            params = c("r_pp", "n"),
            derived_from = "resource parameters")
    )
}

#' Signal the changes to species parameters that cannot take effect
#'
#' Goes through the rate arrays that can be frozen, see [frozen_rate_params()],
#' and raises a [signal_frozen()] condition for each frozen array that one of
#' the changed species parameters feeds. This is what turns "the model no
#' longer follows the species parameters" into a warning the user sees at the
#' moment they make the change. It is one of the diagnostics that only
#' [given_species_params<-()] gives, see there.
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
        # "species parameters" -> "species parameter" for a single one
        plural <- if (is.null(entry$derived_from)) {
            "species parameters"
        } else {
            entry$derived_from
        }
        singular <- sub("s$", "", plural)
        signal_frozen(var, paste0(
            "Your change to the ",
            ngettext(length(affected), singular, plural), " ",
            paste0("`", affected, "`", collapse = ", "),
            " has not taken effect because the ", entry$quantity,
            " has been set manually and so is no longer calculated from the ",
            plural, ". Call `", entry$reset_call, "` if you want the ",
            entry$quantity, " to be calculated from the ", plural, " again."))
    }
    invisible(NULL)
}

# The species parameters that mizer only uses to calculate a default for
# another one, and so ignores once that other one has been given. Each entry is
# named after the parameter that is ignored and gives the parameter(s) that take
# precedence over it.
overridden_species_params <- function() {
    list(f0 = "gamma", fc = "ks", age_mat = "h", k_vb = c("h", "age_mat"))
}

#' Signal the changes that are ignored because another parameter was given
#'
#' Some species parameters are only used to calculate a default for another
#' one: `f0` for `gamma`, `fc` for `ks`, `age_mat` for `h`, and `k_vb` for
#' `h` or `age_mat`. Once the other one has been given, the model no longer
#' consults them, so changing them has no effect. This raises a warning about
#' that. It is one of the diagnostics that only [given_species_params<-()]
#' gives, see there.
#'
#' Only a value that is there can be ignored, so this is asked about the
#' species that were *given* a value, not about every species whose value
#' changed: clearing a value to `NA` is a change, but not one this has anything
#' to say about.
#'
#' @param given The given species parameters, as they are before the change.
#' @param changed A named list with one logical vector per changed column,
#'   saying which species were given a value, as built by
#'   [given_species_params<-()].
#'
#' @return `NULL` invisibly. Called for its side effect of signalling.
#' @concept helper
signal_ignored_changes <- function(given, changed) {
    for (par in names(overridden_species_params())) {
        if (!(par %in% names(changed))) {
            next
        }
        takes_precedence_list <- overridden_species_params()[[par]]
        active_changes <- changed[[par]]
        for (takes_precedence in takes_precedence_list) {
            if (!(takes_precedence %in% names(given))) {
                next
            }
            overridden <- active_changes & !is.na(given[[takes_precedence]])
            if (any(overridden)) {
                signal_info(par, paste0(
                    "The values you specified for `", par, "` will not lead ",
                    "to a re-calculation of `", takes_precedence,
                    "` because its value was explicitly given."),
                    level = 1, severity = "warning", unhandled = "show")
                active_changes[overridden] <- FALSE
            }
        }
    }
    invisible(NULL)
}

#' Signal a gear parameter changed through the given species parameters
#'
#' Mizer looks for the gear parameters in the gear parameter table, which is
#' read only when the model is built, so changing one of them through the
#' species parameters does not reach the model. This is one of the diagnostics
#' that only [given_species_params<-()] gives; [species_params<-()] stays
#' quiet, see there.
#'
#' `yield_observed` is the exception. It belongs in [gear_params()], which gives
#' the observed yield per gear and species, and that is where it should be set,
#' which is what this reports. But it feeds no rate, so a value in the species
#' parameters is not lost: [get_yield_observed()] falls back to it for any
#' species that has no observation among the gear parameters.
#'
#' @param changed A named list with one entry per changed column, or a
#'   character vector of the changed column names.
#'
#' @return `NULL` invisibly. Called for its side effect of signalling.
#' @concept helper
signal_gear_params_changes <- function(changed) {
    changed <- if (is.character(changed)) changed else names(changed)
    gear_pars <- c("catchability", "selectivity", "l50", "l25", "sel_func")
    if (any(gear_pars %in% changed)) {
        signal_info("gear_params", paste0(
            "To make changes to gears you should use `gear_params()<-`, not ",
            "`species_params()`."),
            level = 1, severity = "warning", unhandled = "show")
    }
    if ("yield_observed" %in% changed) {
        signal_info("yield_observed", paste0(
            "To change the observed yield you should use `gear_params()<-`, ",
            "not `species_params()`."),
            level = 1, severity = "warning", unhandled = "show")
    }
    invisible(NULL)
}
