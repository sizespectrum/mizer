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
#' @param fallback The level to use when the option is not set. Defaults to 3,
#'   which reports everything.
#'
#' @return A single number, or `NA` to leave the reporting to a handler further
#'   out.
#' @concept helper
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
#' @param expr The expression to evaluate. It is evaluated in the calling
#'   environment, so assignments made in it have the same effect as they would
#'   have without this wrapper.
#' @param info_level The level of information to report, or `NA` to leave the
#'   reporting to a handler further out. Defaults to
#'   [default_info_level()], which consults the `mizer_info_level` option.
#'
#' @return The value of `expr`.
#' @concept helper
with_info_level <- function(expr, info_level = default_info_level()) {
    # A handler further out is already collecting, or we were told to leave the
    # reporting to one: evaluate without collecting or muffling anything.
    if (info_reporting$active ||
        !is.numeric(info_level) || length(info_level) != 1 ||
        is.na(info_level)) {
        return(expr)
    }
    reports <- list()
    collect_info <- function(cnd) {
        # The fields are defaulted rather than required, so that a condition
        # raised by an extension package that does not use `signal_info()` is
        # still reported rather than throwing.
        level <- if (is.null(cnd$level)) 3 else cnd$level
        severity <- if (is.null(cnd$severity)) "info" else cnd$severity
        if (level <= info_level) {
            key <- paste(severity, cnd$var, cnd$message, sep = "\r")
            reports[[key]] <<- list(severity = severity,
                                    message = cnd$message)
        }
        # Muffle even the conditions that we do not report, because we have
        # taken responsibility for them: `info_level = 0` means silence.
        cnd_muffle(cnd)
    }
    info_reporting$active <- TRUE
    on.exit(info_reporting$active <- FALSE, add = TRUE)
    result <- withCallingHandlers(expr, info_about_default = collect_info)

    severities <- vapply(reports, `[[`, character(1), "severity")
    messages <- vapply(reports, `[[`, character(1), "message")
    if (any(severities == "info")) {
        message(paste(messages[severities == "info"], collapse = "\n"))
    }
    if (any(severities == "warning")) {
        warning(paste(messages[severities == "warning"], collapse = "\n"),
                call. = FALSE)
    }
    result
}

#' Signal information about a choice mizer made
#'
#' Raises the condition that [with_info_level()] collects. This is the way for
#' mizer, and for anything extending it, to tell the user about a default it
#' filled in, an input it adjusted or an instruction it could not carry out,
#' without deciding on its own how loudly to say it: the handler installed by
#' whichever function the user actually called does that.
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
#'   sense as part of a report about a whole model. `"show"` gives the message
#'   anyway.
#' @param class Further classes to give the condition, for code that wants to
#'   catch a particular kind of report.
#'
#' @return `NULL` invisibly. Called for its side effect of signalling.
#' @concept helper
signal_info <- function(var, message, level = 3,
                        severity = c("info", "warning"),
                        unhandled = c("drop", "show"),
                        class = character()) {
    severity <- match.arg(severity)
    unhandled <- match.arg(unhandled)
    emit <- if (unhandled == "show") inform else signal
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
#' @concept helper
signal_not_recalculated <- function(var, quantity, reset_call,
                                    derived_from = "species parameters") {
    signal_info(var, paste0(
        "The ", quantity, " has been set manually and so is not recalculated ",
        "from the ", derived_from, ". Call `", reset_call, "` if you want the ",
        quantity, " to be calculated from the ", derived_from, " again."),
        level = 1, unhandled = "show")
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
