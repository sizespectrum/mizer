# Cache resource arrays once to avoid recomputing them in every test block.
mort_small <- getResourceMort(NS_params_small)
n_resource_small <- NResource(NS_sim_small)

test_that("ArrayResourceBySize constructor works", {
    vec <- 1:10
    names(vec) <- signif(NS_params_small@w_full[1:10], 3)

    rate <- ArrayResourceBySize(vec, value_name = "Test rate", units = "1/year")

    expect_s3_class(rate, "ArrayResourceBySize")
    expect_true(is.ArrayResourceBySize(rate))
    expect_false(is.ArrayResourceBySize(vec))
    expect_identical(length(rate), length(vec))
    expect_identical(attr(rate, "value_name"), "Test rate")
    expect_identical(attr(rate, "units"), "1/year")
})

test_that("ArrayResourceBySize constructor validates input", {
    expect_error(ArrayResourceBySize(matrix(1:4, 2, 2)), "`x` must be a numeric vector.")
})

test_that("Rate functions return ArrayResourceBySize", {
    expect_true(is.ArrayResourceBySize(mort_small))
    expect_true(is.ArrayResourceBySize(initialNResource(NS_params_small)))
})

test_that("print.ArrayResourceBySize works", {
    expect_output(print(mort_small), "Resource mortality")
    expect_output(print(mort_small), "sizes")
    expect_output(print(mort_small), "1/year")
})

test_that("print.ArrayResourceBySize truncates a wide size grid", {
    mort <- getResourceMort(NS_params)
    expect_output(print(mort), "showing .* of \\d+ sizes")
    expect_output(print(mort), "log-spaced")
    expect_output(print(mort), "as.data.frame")
})

test_that("summary.ArrayResourceBySize works", {
    s <- summary(mort_small)
    expect_s3_class(s, "summary.ArrayResourceBySize")
    expect_identical(s$value_name, "Resource mortality")
    expect_output(print(s), "Resource mortality")
})

test_that("str.ArrayResourceBySize works", {
    expect_output(str(mort_small), "ArrayResourceBySize")
    expect_output(str(mort_small), "Resource mortality")
})

test_that("ArrayTimeByResourceBySize constructor works", {
    mat <- matrix(1:20, nrow = 2, ncol = 10)
    rownames(mat) <- as.character(2000:2001)
    colnames(mat) <- signif(NS_params_small@w_full[1:10], 3)

    out <- ArrayTimeByResourceBySize(mat, value_name = "Test rate", units = "1/g")

    expect_s3_class(out, "ArrayTimeByResourceBySize")
    expect_true(is.ArrayTimeByResourceBySize(out))
    expect_false(is.ArrayTimeByResourceBySize(mat))
    expect_identical(dim(out), dim(mat))
    expect_identical(attr(out, "value_name"), "Test rate")
    expect_identical(attr(out, "units"), "1/g")
})

test_that("ArrayTimeByResourceBySize constructor validates input", {
    expect_error(ArrayTimeByResourceBySize(1:10), "`x` must be a matrix.")
})

test_that("NResource returns ArrayTimeByResourceBySize", {
    expect_true(is.ArrayTimeByResourceBySize(n_resource_small))
})

test_that("print.ArrayTimeByResourceBySize works", {
    expect_output(print(n_resource_small), "Number density")
    expect_output(print(n_resource_small), "times x")
    expect_output(print(n_resource_small), "1/g")
})

test_that("print.ArrayTimeByResourceBySize truncates a long time series", {
    sim_long <- suppressMessages(
        project(NS_params_small, t_max = 60, dt = 0.5, t_save = 0.5,
               progress_bar = FALSE))
    n_resource <- NResource(sim_long)
    expect_identical(nrow(n_resource), 121L)
    out <- paste(capture.output(print(n_resource)), collapse = "\n")
    expect_match(out, "\\b0\\b")
    expect_match(out, "\\b60\\b")
    expect_match(out, "showing 8 of 121 times")
    expect_match(out, "evenly spaced")
})

test_that("summary.ArrayTimeByResourceBySize works", {
    s <- summary(n_resource_small)
    expect_s3_class(s, "summary.ArrayTimeByResourceBySize")
    expect_identical(s$value_name, "Number density")
    expect_output(print(s), "Number density")
})

test_that("str.ArrayTimeByResourceBySize works", {
    expect_output(str(n_resource_small), "ArrayTimeByResourceBySize")
    expect_output(str(n_resource_small), "Number density")
})
