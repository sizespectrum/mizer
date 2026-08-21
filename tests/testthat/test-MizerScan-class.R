# Fast settings, so that the integration tests stay cheap. Note that t_per must
# be a whole multiple of both dt and t_save.
fast <- list(t_max = 6, t_per = 1.5, dt = 0.5, t_save = 0.5, t_sample = 3,
             progress_bar = FALSE)

scan_fast <- function(params, ...) {
    do.call(scanModel, c(list(params = params, ...), fast))
}

# Where the maximum is -------------------------------------------------------

test_that("at_max records where each series is largest, and stays fresh", {
    scan <- suppressMessages(suppressWarnings(
        scanModel(NS_params_steady_small, scan_values = c(0.4, 0.5, 0.6),
                  set_func = scanEffort("Otter"), value_func = getYield,
                  progress_bar = FALSE)))
    at_max <- attr(scan, "at_max")
    max_value <- attr(scan, "max_value")
    expect_identical(names(at_max), unique(as.character(scan$Species)))
    for (s in names(at_max)) {
        rows <- as.character(scan$Species) == s
        expect_equal(max_value[[s]], max(scan[[2]][rows]))
        expect_equal(at_max[[s]], scan[[1]][rows][[which.max(scan[[2]][rows])]])
    }
    # Subsetting recomputes rather than copies, so a subset cannot carry a
    # maximum taken over rows it no longer contains.
    sub <- scan[scan[[1]] <= 0.5, ]
    expect_true(all(attr(sub, "at_max") <= 0.5, na.rm = TRUE))
})

# Methods --------------------------------------------------------------------

test_that("the MizerScan methods behave", {
    scan <- suppressMessages(suppressWarnings(
        scan_fast(NS_params_small, scan_values = c(0, 1),
                  set_func = scanEffort())))

    expect_output(print(scan), "Biomass")
    expect_output(print(scan), "Fishing effort")

    s <- summary(scan)
    expect_s3_class(s, "summary.MizerScan")
    expect_identical(nrow(s$per_series), 3L)
    expect_output(print(s), "at_max")

    # str() must not dump the whole MizerParams
    out <- capture.output(str(scan))
    expect_true(any(grepl("MizerScan", out)))
    expect_false(any(grepl("intake_max", out)))

    # `[` keeps the class while the result is still a scan, and drops it once
    # the result no longer has the three variables a scan needs.
    expect_true(is.MizerScan(scan[scan$Species == "Cod", ]))
    expect_false(is.MizerScan(scan[, 1:2]))
    expect_type(scan[["ymin"]], "double")

    plain <- as.data.frame(scan)
    expect_false(is.MizerScan(plain))
    expect_null(attr(plain, "params"))
    expect_null(attr(plain, "value_name"))
})

test_that("MizerScan() checks its input", {
    expect_error(MizerScan("a"), "must be a data frame")
    expect_error(MizerScan(data.frame(a = 1, b = 2)), "at least three columns")
    expect_error(MizerScan(data.frame(a = 1, b = 2, c = 3)),
                 "`Species` column")
})

test_that("plot.MizerScan() draws all three styles", {
    scan <- suppressMessages(suppressWarnings(
        scan_fast(NS_params_small, scan_values = c(0, 0.5, 1),
                  set_func = scanEffort())))
    for (style in c("ribbon", "envelope", "line")) {
        p <- plot(scan, style = style)
        expect_s3_class(p, "mizer_plot")
        expect_s3_class(p, "ggplot")
    }
    expect_s3_class(plot(scan, return_data = TRUE), "data.frame")
    # Reference lines can be given explicitly as well as stored
    expect_s3_class(plot(scan, reference_lines = c(F_MSY = 0.5)), "ggplot")
    expect_s3_class(plot(scan, mark_max = TRUE), "ggplot")

    # A scan that has lost its params, as dplyr::filter() would leave it,
    # cannot be plotted and says why.
    stripped <- scan
    attr(stripped, "params") <- NULL
    expect_error(plot(stripped), "lost its `params` attribute")
})
