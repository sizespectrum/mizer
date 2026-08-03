test_that("pick_log_spaced_indices returns everything below the threshold", {
    expect_identical(pick_log_spaced_indices(5, max_n = 10), 1:5)
    # `threshold` delays truncation beyond max_n alone
    expect_identical(pick_log_spaced_indices(15, max_n = 10, threshold = 20),
                     1:15)
})

test_that("pick_log_spaced_indices thins to at most max_n sorted unique indices", {
    idx <- pick_log_spaced_indices(100, max_n = 10)
    expect_length(idx, 10)
    expect_identical(idx, sort(unique(idx)))
    # the ends of the range are always kept
    expect_identical(range(idx), c(1, 100))
})

test_that("pick_head_indices takes a plain leading subset", {
    expect_identical(pick_head_indices(5, head_n = 10), 1:5)
    expect_identical(pick_head_indices(100, head_n = 3), 1:3)
    expect_identical(pick_head_indices(15, head_n = 3, threshold = 20), 1:15)
})

test_that("matrix_display_width counts row labels and all columns", {
    mat <- matrix(1:4, nrow = 2, dimnames = list(c("aa", "b"), NULL))
    formatted <- format(mat)
    expected <- max(nchar(rownames(mat))) +
        sum(apply(formatted, 2, function(col) max(nchar(col))) + 1)
    expect_identical(matrix_display_width(mat), expected)
})

test_that("vector_display_width uses the wider of name and value", {
    vec <- c(alpha = 1, b = 22)
    # "alpha" (5) is wider than its formatted value, "22" ties with its name
    expect_identical(vector_display_width(vec),
                     sum(pmax(nchar(names(vec)), nchar(format(vec))) + 1))
})

test_that("fit_log_spaced_k returns the largest k that fits the width", {
    # width_fn(k) = 10 * k, so with width 55 the largest fitting k is 5
    expect_identical(fit_log_spaced_k(100, max_k = 10, min_k = 1,
                                      width_fn = function(k) 10 * k,
                                      width = 55), 5L)
})

test_that("fit_log_spaced_k is bounded by n and by min_k", {
    expect_identical(fit_log_spaced_k(3, max_k = 10, min_k = 1,
                                      width_fn = function(k) 1,
                                      width = 80), 3L)
    # nothing fits, so it falls back to min_k
    expect_identical(fit_log_spaced_k(100, max_k = 10, min_k = 2,
                                      width_fn = function(k) 1000,
                                      width = 80), 2)
})

test_that("parse_numeric_labels parses numeric labels and falls back otherwise", {
    expect_identical(parse_numeric_labels(c("1", "2.5", "10"), 3), c(1, 2.5, 10))
    expect_identical(parse_numeric_labels(NULL, 3), 1:3)
    expect_identical(parse_numeric_labels(c("a", "b"), 2), 1:2)
})

test_that("format_truncation_note reports the counts and optional detail", {
    expect_identical(format_truncation_note(10, 226, "sizes"),
                     "... showing 10 of 226 sizes; use as.data.frame() for the full data.")
    expect_identical(format_truncation_note(10, 226, "sizes", "0.001-4e+04 g"),
                     "... showing 10 of 226 sizes (0.001-4e+04 g); use as.data.frame() for the full data.")
    # an empty detail is treated like no detail at all
    expect_identical(format_truncation_note(1, 2, "times", ""),
                     "... showing 1 of 2 times; use as.data.frame() for the full data.")
})

test_that("format_size_range_detail reports the weight range in grams", {
    expect_identical(format_size_range_detail(c(0.0012345, 1, 40000)),
                     "0.00123-40000 g, log-spaced")
})

test_that("format_time_range_detail reports the time range", {
    expect_identical(format_time_range_detail(c(1, 5, 10)),
                     "1-10, evenly spaced")
})
