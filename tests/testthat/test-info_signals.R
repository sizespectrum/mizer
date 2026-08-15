# with_info_level() ----

test_that("with_info_level() collects information into a single message", {
    expect_message(
        with_info_level({
            signal_info("a", "first", level = 1)
            signal_info("b", "second")
        }),
        "first\nsecond")
})

test_that("with_info_level() collapses repeats but keeps distinct reports", {
    # The same report twice takes up one line ...
    expect_message(
        with_info_level({
            signal_info("a", "same", level = 1)
            signal_info("a", "same", level = 1)
        }),
        "^same\n$")
    # ... but two things said about one quantity are both kept.
    expect_message(
        with_info_level({
            signal_info("a", "first", level = 1)
            signal_info("a", "second", level = 1)
        }),
        "^first\nsecond\n$")
})

test_that("with_info_level() drops information above the info level", {
    expect_message(
        with_info_level(info_level = 1, {
            signal_info("a", "important", level = 1)
            signal_info("b", "chatter", level = 3)
        }),
        "^important\n$")
    # info_level = 0 is silence
    expect_silent(
        with_info_level(info_level = 0, signal_info("a", "important",
                                                    level = 1)))
    # even for a signal that shows itself when unhandled
    expect_silent(
        with_info_level(info_level = 0,
                        signal_info("a", "important", level = 1,
                                    unhandled = "show")))
})

test_that("with_info_level() evaluates the expression in the calling frame", {
    result <- suppressMessages(with_info_level({
        x <- 2
        signal_info("a", "info", level = 1)
        x + 1
    }))
    expect_identical(result, 3)
    expect_identical(x, 2)
})

test_that("with_info_level() handlers nest, the outermost one reporting", {
    # The inner handler steps aside, so there is one report, not two.
    expect_message(
        with_info_level({
            with_info_level(signal_info("a", "inner", level = 1))
        }),
        "^inner\n$")
    # The inner info level does not apply, the outer one does
    expect_silent(
        with_info_level(info_level = 0, {
            with_info_level(info_level = 3, signal_info("a", "inner"))
        }))
    # but `info_level = 0` silences its own expression even so
    expect_silent(
        with_info_level({
            with_info_level(info_level = 0, signal_info("a", "inner",
                                                        level = 1))
        }))
    # and the outer handler carries on reporting afterwards
    expect_message(
        with_info_level({
            with_info_level(info_level = 0, signal_info("a", "quiet",
                                                        level = 1))
            signal_info("b", "loud", level = 1)
        }),
        "^loud\n$")
    # `NA` asks for the same deferral explicitly
    expect_message(
        with_info_level({
            with_info_level(info_level = NA, signal_info("a", "inner",
                                                         level = 1))
        }),
        "^inner\n$")
})

test_that("with_info_level() reports when the expression returns early", {
    f <- function() {
        with_info_level({
            signal_info("a", "early", level = 1)
            return("returned")
            signal_info("b", "never", level = 1)  # nocov
        })
    }
    expect_message(expect_identical(f(), "returned"), "^early\n$")
})

test_that("with_info_level() releases the nesting flag when expr fails", {
    expect_error(with_info_level(stop("boom")), "boom")
    expect_false(info_reporting$active)
    expect_message(with_info_level(signal_info("a", "after", level = 1)),
                   "^after\n$")
})

test_that("the mizer_info_level option sets the default", {
    withr::local_options(mizer_info_level = 0)
    expect_identical(default_info_level(), 0)
    expect_identical(default_info_level(2), 0)
    expect_silent(with_info_level(signal_info("a", "chatter", level = 1)))
    # and it reaches the functions that have no info_level argument
    params <- NS_params_small
    metab(params) <- metab(params)
    sp <- species_params(params)
    sp$ks <- sp$ks * 2
    expect_silent(species_params(params) <- sp)
})

test_that("default_info_level() falls back when the option is not set", {
    withr::local_options(mizer_info_level = NULL)
    expect_identical(default_info_level(), 3)
    expect_identical(default_info_level(2), 2)
})

# signal_info() ----

test_that("signal_info() says nothing when unhandled unless asked to", {
    expect_silent(signal_info("a", "chatter"))
    expect_message(signal_info("a", "chatter", unhandled = "show"),
                   "^chatter$")
    # An unhandled report keeps its severity
    expect_warning(signal_info("a", "alarm", severity = "warning",
                               unhandled = "show"),
                   "^alarm$")
    expect_silent(signal_info("a", "alarm", severity = "warning"))
})

test_that("signal_info() reports at the requested severity", {
    expect_message(with_info_level(signal_info("a", "note", level = 1)),
                   "^note\n$")
    expect_warning(
        with_info_level(signal_info("a", "alarm", level = 1,
                                    severity = "warning")),
        "^alarm$")
})

test_that("signal_info() rejects an unknown severity", {
    expect_error(signal_info("a", "x", severity = "fatal"), "'arg'")
})

test_that("with_info_level() copes with a condition that has no fields", {
    # An extension package may raise the condition itself
    expect_message(
        with_info_level(signal("bare", class = "info_about_default")),
        "^bare\n$")
    expect_silent(
        with_info_level(info_level = 0,
                        signal("bare", class = "info_about_default")))
})

# signal_frozen() ----

test_that("with_info_level() reports frozen signals as a warning", {
    expect_warning(
        with_info_level(signal_frozen("metab", "frozen")),
        "^frozen$")
    # It is a warning and not a message, so that it survives
    # `suppressMessages()`, see `species_params<-()`.
    expect_warning(
        suppressMessages(with_info_level(signal_frozen("metab", "frozen"))),
        "^frozen$")
    # Messages and warnings are reported separately
    expect_warning(
        expect_message(
            with_info_level({
                signal_info("a", "info", level = 1)
                signal_frozen("metab", "frozen")
            }),
            "^info\n$"),
        "^frozen$")
    # and `info_level = 0` still silences both
    expect_silent(
        with_info_level(info_level = 0, signal_frozen("metab", "frozen")))
})

test_that("signal_frozen() surfaces as a warning when nobody is listening", {
    expect_warning(signal_frozen("metab", "frozen"), "^frozen$")
})

test_that("signal_not_recalculated() builds a message naming the way back", {
    # Surfaces as a message when nobody is listening
    expect_message(
        signal_not_recalculated("metab", "metabolic rate",
                                "setMetabolicRate(params, reset = TRUE)"),
        "The metabolic rate has been set manually")
    expect_message(
        with_info_level(
            signal_not_recalculated("metab", "metabolic rate",
                                    "setMetabolicRate(params, reset = TRUE)")),
        "The metabolic rate has been set manually.*not recalculated from the species parameters.*setMetabolicRate\\(params, reset = TRUE\\)")
    expect_message(
        with_info_level(
            signal_not_recalculated("selectivity", "selectivity",
                                    "setFishing(params, reset = TRUE)",
                                    derived_from = "gear parameters")),
        "not recalculated from the gear parameters")
})

# signal_frozen_changes() ----

test_that("signal_frozen_changes() only signals about frozen quantities", {
    params <- NS_params_small
    expect_silent(signal_frozen_changes(params, c("ks", "gamma", "w_mat")))

    metab(params) <- metab(params)
    # `ks` feeds the metabolic rate, `beta` does not
    expect_warning(with_info_level(signal_frozen_changes(params, "ks")),
                   "Your change to the species parameter `ks`.*metabolic rate")
    expect_silent(with_info_level(signal_frozen_changes(params, "beta")))
    # and nothing at all is signalled when nothing changed
    expect_silent(with_info_level(signal_frozen_changes(params, character(0))))
})

test_that("signal_frozen_changes() lists all the affected parameters", {
    params <- NS_params_small
    metab(params) <- metab(params)
    expect_warning(
        with_info_level(signal_frozen_changes(params, c("ks", "p", "beta"))),
        "species parameters `ks`, `p` has not taken effect")
})

# signal_ignored_changes() ----

test_that("signal_ignored_changes() warns about a parameter that is overruled", {
    given <- given_species_params(NS_params_small)
    all_sp <- rep(TRUE, nrow(given))
    # `gamma` is given, so a change to `f0` has no effect
    expect_true(all(!is.na(given$gamma)))
    expect_warning(
        with_info_level(signal_ignored_changes(given, list(f0 = all_sp))),
        "values for `f0` that are going to be ignored because values for `gamma`")
    # but not for a species whose `gamma` is not given
    given$gamma[[1]] <- NA
    expect_silent(
        with_info_level(signal_ignored_changes(given, list(f0 = c(TRUE, FALSE, FALSE)))))
    expect_warning(
        with_info_level(signal_ignored_changes(given, list(f0 = c(FALSE, TRUE, FALSE)))),
        "`f0`")
    # and not at all when the overruling parameter was never given
    given$gamma <- NULL
    expect_silent(
        with_info_level(signal_ignored_changes(given, list(f0 = all_sp))))
})

test_that("signal_gear_params_changes() warns about gear parameters", {
    expect_warning(
        with_info_level(signal_gear_params_changes(list(l50 = TRUE))),
        "you should use `gear_params\\(\\)<-`")
    expect_warning(
        with_info_level(signal_gear_params_changes("yield_observed")),
        "observed yield")
    expect_silent(with_info_level(signal_gear_params_changes("gamma")))
})

test_that("only the given species params report a gear parameter change", {
    quiet <- NS_params_small
    expect_silent(species_params(quiet)$catchability <- 2)
    loud <- NS_params_small
    expect_warning(given_species_params(loud)$catchability <- 2,
                   "you should use `gear_params\\(\\)<-`")
})

test_that("only the given species params report a `yield_observed` change", {
    # `yield_observed` belongs in `gear_params()`. Only the setter that carries
    # the diagnostics says so; `species_params<-()` is the quiet one.
    quiet <- NS_params_small
    expect_silent(species_params(quiet)$yield_observed <- c(1, 2, 3))
    loud <- NS_params_small
    expect_warning(given_species_params(loud)$yield_observed <- c(4, 5, 6),
                   "observed yield")
})
