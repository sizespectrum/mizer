# with_info_level() ----

test_that("with_info_level() collects information into a single message", {
    expect_message(
        with_info_level({
            signal("first", class = "info_about_default", var = "a", level = 1)
            signal("second", class = "info_about_default", var = "b", level = 3)
        }),
        "first\nsecond")
})

test_that("with_info_level() reports only one message per var", {
    expect_message(
        with_info_level({
            signal("old", class = "info_about_default", var = "a", level = 1)
            signal("new", class = "info_about_default", var = "a", level = 1)
        }),
        "^new\n$")
})

test_that("with_info_level() drops information above the info level", {
    expect_message(
        with_info_level(info_level = 1, {
            signal("important", class = "info_about_default", var = "a",
                   level = 1)
            signal("chatter", class = "info_about_default", var = "b",
                   level = 3)
        }),
        "^important\n$")
    # info_level = 0 is silence
    expect_silent(
        with_info_level(info_level = 0, {
            signal("important", class = "info_about_default", var = "a",
                   level = 1)
        }))
})

test_that("with_info_level() evaluates the expression in the calling frame", {
    result <- suppressMessages(with_info_level({
        x <- 2
        signal("info", class = "info_about_default", var = "a", level = 1)
        x + 1
    }))
    expect_identical(result, 3)
    expect_identical(x, 2)
})

test_that("with_info_level(info_level = NA) defers to the outer handler", {
    # The inner call neither reports nor swallows, so the outer one reports.
    expect_message(
        with_info_level({
            with_info_level(info_level = NA, {
                signal("inner", class = "info_about_default", var = "a",
                       level = 1)
            })
        }),
        "^inner\n$")
})

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
                signal("info", class = "info_about_default", var = "a",
                       level = 1)
                signal_frozen("metab", "frozen")
            }),
            "^info\n$"),
        "^frozen$")
    # and `info_level = 0` still silences both
    expect_silent(
        with_info_level(info_level = 0, signal_frozen("metab", "frozen")))
})

# signal_frozen() ----

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
