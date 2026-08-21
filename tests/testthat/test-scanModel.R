# Pure helpers ---------------------------------------------------------------

test_that("sweep_arms() works outwards from the current value", {
    v <- c(0, 0.5, 1, 1.5)
    # Without a current value there is one arm, in the order given, which is
    # what lets a decreasing vector trace a hysteresis branch deliberately.
    expect_identical(mizer:::sweep_arms(v), list(seq_along(v)))
    expect_identical(mizer:::sweep_arms(rev(v)), list(seq_along(v)))
    # From 1 there are two arms: down through the smaller values, and up.
    arms <- mizer:::sweep_arms(v, 1)
    expect_length(arms, 2)
    expect_identical(v[arms[[1]]], c(0.5, 0))
    expect_identical(v[arms[[2]]], c(1, 1.5))
    # A value sitting exactly at the current one goes into the ascending arm,
    # and an empty arm is dropped rather than returned.
    expect_identical(mizer:::sweep_arms(v, 0), list(seq_along(v)))
    expect_identical(lapply(mizer:::sweep_arms(v, 2), function(a) v[a]),
                     list(c(1.5, 1, 0.5, 0)))
})

test_that("each arm of the scan restarts from the model as given", {
    # The arms must not be run as one sequence: doing so would start the
    # ascending arm at the attractor reached at the far end of the descending
    # arm, which in a hysteretic model follows the wrong branch.
    seen <- list()
    spy <- function(params, value) {
        seen[[length(seen) + 1L]] <<- list(value = value,
                                           n = params@initial_n)
        scanEffort()(params, value)
    }
    v <- c(0, 0.5, 1, 1.5)
    suppressMessages(suppressWarnings(
        scan_fast(NS_params_small, scan_values = v, set_func = spy,
                  current_scan_value = 1)))

    values <- vapply(seen, function(x) x$value, numeric(1))
    expect_identical(values, c(0.5, 0, 1, 1.5))
    # The first call of each arm sees the unmodified model ...
    original <- NS_params_small@initial_n
    expect_equal(seen[[1]]$n, original)   # start of the descending arm
    expect_equal(seen[[3]]$n, original)   # start of the ascending arm
    # ... while continuation within an arm does not.
    expect_false(isTRUE(all.equal(seen[[2]]$n, original)))
})

test_that("as_series_matrix() normalises what a value function returns", {
    m <- matrix(1:6, nrow = 3, dimnames = list(NULL, c("a", "b")))
    expect_identical(mizer:::as_series_matrix(m), m)

    v <- mizer:::as_series_matrix(c(`0` = 1, `1` = 2), default_name = "Q")
    expect_true(is.matrix(v))
    expect_identical(dim(v), c(2L, 1L))
    expect_identical(colnames(v), "Q")

    unnamed <- mizer:::as_series_matrix(matrix(1:4, nrow = 2))
    expect_identical(colnames(unnamed), c("Series 1", "Series 2"))

    expect_error(mizer:::as_series_matrix(data.frame(slope = 1)),
                 "wrap the column|pick out the column")
    expect_error(mizer:::as_series_matrix(array(1, dim = c(2, 2, 2))),
                 "3 dimensions")
    expect_error(mizer:::as_series_matrix("a"), "must return numbers")
})

test_that("the last sample of a whole period is left out of the mean", {
    # A signal sampled over exactly one period has a final sample that repeats
    # the phase of the first, so including it biases the mean.
    # A cosine, not a sine: a sine sampled over a whole period is symmetric
    # about zero, so the duplicated endpoint would not bias its mean and the
    # test would pass whatever the code did.
    n <- 20
    y <- cos(2 * pi * (0:n) / n)
    expect_equal(mean(y[seq_len(n)]), 0)
    expect_gt(abs(mean(y)), 1e-3)
})
test_that("scanModel() insists that set_func returns a MizerParams", {
    expect_error(scan_fast(NS_params_small, scan_values = c(1, 2),
                           set_func = function(params, value) 42),
                 "must return a MizerParams")
})

test_that("params_as_sim() carries the whole state, including the effort", {
    sim <- mizer:::params_as_sim(NS_params_small)
    # getYield() reads sim@effort, which MizerSim() would leave as NA.
    expect_equal(unclass(getYield(sim))[1, ], unclass(getYield(NS_params_small)),
                 ignore_attr = TRUE)
    expect_equal(unclass(getBiomass(sim))[1, ],
                 unclass(getBiomass(NS_params_small)), ignore_attr = TRUE)
    expect_true(all(is.finite(getSSB(sim))))
})

# The scan itself ------------------------------------------------------------

test_that("scanModel() returns a MizerScan laid out as documented", {
    scan <- suppressMessages(suppressWarnings(
        scan_fast(NS_params_small, scan_values = c(0, 0.5, 1),
                  set_func = scanEffort())))
    expect_true(is.MizerScan(scan))
    expect_s3_class(scan, "data.frame")
    expect_identical(names(scan),
                     c("Fishing effort", "Biomass", "Species", "ymin", "ymax",
                       "type", "settled", "period", "residual"))
    expect_identical(nrow(scan), 9L)
    # The metadata is lifted off what the value function returned.
    expect_identical(attr(scan, "value_name"), "Biomass")
    expect_identical(attr(scan, "value_units"), "g")
    expect_identical(attr(scan, "scan_name"), "Fishing effort")
    expect_true(all(scan$ymax >= scan$ymin))
    expect_true(all(scan[[2]] >= scan$ymin & scan[[2]] <= scan$ymax))
    # Rows come back in the order the scan values were given.
    expect_identical(unique(scan[[1]]), c(0, 0.5, 1))
})

test_that("a fixed point is not sampled at all", {
    scan <- suppressMessages(suppressWarnings(
        scanModel(NS_params_steady_small, scan_values = c(0.4, 0.5),
                  set_func = scanEffort("Otter"), value_func = getYield,
                  progress_bar = FALSE)))
    fixed <- scan$type == "below_tolerance"
    expect_true(any(fixed))
    # Identical, not merely equal: the value is read off a single snapshot.
    expect_identical(scan[[2]][fixed], scan$ymin[fixed])
    expect_identical(scan[[2]][fixed], scan$ymax[fixed])
    expect_true(all(is.na(scan$period[fixed])))
})

test_that("scanModel() accepts a value function returning a plain vector", {
    scan <- suppressMessages(suppressWarnings(
        scan_fast(NS_params_small, scan_values = c(0, 1),
                  set_func = scanEffort(),
                  value_func = function(sim) getMeanWeight(sim),
                  value_name = "Mean weight", value_units = "g")))
    expect_identical(unique(as.character(scan$Species)), "Mean weight")
    expect_identical(attr(scan, "value_name"), "Mean weight")
    # A series that is not a species must still be drawn, not silently dropped
    # by plotDataFrame() for want of a colour.
    p <- plot(scan)
    expect_gt(nrow(ggplot2::layer_data(p, 1)), 0)
})

test_that("series that are not species can be selected", {
    scan <- suppressMessages(suppressWarnings(
        scan_fast(NS_params_small, scan_values = c(0, 1),
                  set_func = scanEffort(),
                  value_func = function(sim) getMeanWeight(sim),
                  value_name = "Mean weight")))
    # Selecting by name must go through the scan's own series, not through
    # valid_species_arg(), which would reject a series that is not a species.
    p <- plot(scan, species = "Mean weight")
    expect_s3_class(p, "ggplot")
    d <- plot(scan, species = "Mean weight", return_data = TRUE)
    expect_identical(unique(as.character(d$Species)), "Mean weight")
    # By index into the series too
    expect_identical(
        unique(as.character(plot(scan, species = 1, return_data = TRUE)$Species)),
        "Mean weight")
    expect_error(plot(scan, species = "Cod"), "no series called")
    expect_error(plot(scan, species = 5), "out of range")
    # A fractional index would otherwise be truncated silently, and an empty
    # selection would quietly return an empty plot.
    expect_error(plot(scan, species = 1.5), "whole numbers")
    expect_error(plot(scan, species = NA_real_), "must not contain NA")
    expect_error(plot(scan, species = NA_character_), "must not contain NA")
    expect_error(plot(scan, species = numeric(0)), "selects no series")
    expect_error(plot(scan, species = character(0)), "selects no series")
    expect_error(plot(scan, species = Inf), "whole numbers")
})

test_that("scanModel() can restrict the scan to some species", {
    scan <- suppressMessages(suppressWarnings(
        scan_fast(NS_params_small, scan_values = c(0, 1),
                  set_func = scanEffort(), species = "Cod")))
    expect_identical(unique(as.character(scan$Species)), "Cod")
    expect_error(suppressMessages(suppressWarnings(
        scan_fast(NS_params_small, scan_values = c(0, 1),
                  set_func = scanEffort(), species = "Nonexistent"))),
        "no series called")
})

test_that("the sweep order does not change the result", {
    args <- list(NS_params_small, scan_values = c(0, 0.5, 1),
                 set_func = scanEffort())
    a <- suppressMessages(suppressWarnings(do.call(scan_fast, args)))
    b <- suppressMessages(suppressWarnings(
        do.call(scan_fast, c(args, list(current_scan_value = 0.5)))))
    expect_identical(a[[1]], b[[1]])
    expect_identical(a$Species, b$Species)
})

test_that("scanModel() reports the scan values that did not settle", {
    expect_message(
        suppressWarnings(scan_fast(NS_params_small, scan_values = c(0, 1),
                                   set_func = scanEffort())),
        "did not settle onto an attractor within 6 years")
})

test_that("current_scan_value = 'auto' asks the setter", {
    scan <- suppressMessages(suppressWarnings(
        scan_fast(NS_params_small, scan_values = c(0.4, 0.6),
                  set_func = scanEffort("Otter"),
                  current_scan_value = "auto")))
    expect_true(is.MizerScan(scan))
    expect_error(suppressMessages(suppressWarnings(
        scan_fast(NS_params_small, scan_values = c(0, 1),
                  set_func = function(params, value) params,
                  current_scan_value = "auto"))),
        "says where the model currently sits")
})

test_that("current_scan_value says so when it cannot matter", {
    # Raised through signal_info() at severity "warning", so that it survives
    # the suppressMessages() the setter paths run.
    expect_warning(suppressMessages(
        scan_fast(NS_params_small, scan_values = c(0, 1),
                  set_func = scanEffort(), continuation = FALSE,
                  current_scan_value = 0.5)),
        "It will be ignored")
})

test_that("a mistyped species is caught before the whole scan has run", {
    # The series are known after the first scan value, so the error must not
    # wait until every projection has been paid for.
    seen <- 0L
    spy <- function(params, value) {
        seen <<- seen + 1L
        scanEffort()(params, value)
    }
    expect_error(suppressMessages(suppressWarnings(
        scan_fast(NS_params_small, scan_values = c(0, 0.5, 1, 1.5),
                  set_func = spy, species = "Nonexistent"))),
        "no series called")
    expect_identical(seen, 1L)
})

test_that("the generic and its method keep the same signature", {
    # The two are written out separately, as mizer's S3 generics are, so that
    # the man page shows the defaults. Nothing else keeps them in step: the
    # generic's defaults are never evaluated, so a drifted one would show up
    # only as a wrong default in the documentation.
    expect_identical(names(formals(scanModel)),
                     names(formals(mizer:::scanModel.MizerParams)))
    expect_identical(
        vapply(formals(scanModel), deparse1, character(1)),
        vapply(formals(mizer:::scanModel.MizerParams), deparse1, character(1)))
})

# Slow: the correctness claim ------------------------------------------------

test_that("a limit cycle is averaged over exactly one period", {
    skip_on_cran()
    skip_unless_experimental()

    # The limit cycle fixture used in test-steady.R
    p_cycle <- suppressMessages(
        projectToSteady(NS_params, effort = 2, t_max = 200, t_per = 1.5,
                        dt = 0.1, tol = 1e-8, info_level = 0))
    expect_identical(attr(p_cycle, "convergence")$type, "cycle")

    scan <- suppressMessages(suppressWarnings(
        scanModel(p_cycle, scan_values = 1,
                  set_func = scanFishingMortality("Herring"),
                  value_func = getYield, species = "Herring",
                  progress_bar = FALSE)))
    expect_identical(scan$type, "cycle")
    expect_gt(scan$period, 0)
    # The yield oscillates, so the average lies strictly inside the range.
    expect_lt(scan$ymin, scan[[2]])
    expect_lt(scan[[2]], scan$ymax)

    # The average over one period agrees with an average over twenty.
    settled <- suppressMessages(suppressWarnings(
        projectToSteady(scanFishingMortality("Herring")(p_cycle, 1),
                        t_max = 100, tol = 0.001, info_level = 0)))
    period <- attr(settled, "convergence")$period
    long <- project(settled, t_max = round(20 * period / 0.1) * 0.1,
                    dt = 0.1, t_save = 0.1, progress_bar = FALSE)
    y <- as.vector(getYield(long)[, "Herring"])
    expect_equal(scan[[2]], mean(y[-length(y)]), tolerance = 0.02)
})
