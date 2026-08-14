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
    # `NA` asks for the same deferral explicitly
    expect_message(
        with_info_level({
            with_info_level(info_level = NA, signal_info("a", "inner",
                                                         level = 1))
        }),
        "^inner\n$")
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

test_that("signal_frozen() surfaces as a message when nobody is listening", {
    expect_message(signal_frozen("metab", "frozen"), "^frozen$")
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

# changed_species_params() ----

test_that("changed_species_params() finds the changed columns", {
    sp <- species_params(NS_params_small)
    new <- sp
    new$ks <- new$ks * 2
    new$new_col <- 1
    expect_setequal(changed_species_params(new, sp), c("ks", "new_col"))
    expect_identical(changed_species_params(sp, sp), character(0))
})
