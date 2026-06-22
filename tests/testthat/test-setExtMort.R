params <- NS_params_small

## setExtMort ----
test_that("setExtMort works", {
    params <- NS_params_small
    params@species_params$z0 <- 2 * params@species_params$z0
    p2 <- setExtMort(params)
    expect_equal(p2@mu_b, 2 * params@mu_b, ignore_attr = TRUE)
    # only mu_b changed
    p2_compare <- p2
    p2_compare@mu_b <- params@mu_b
    p2_compare@species_params[c("d", "z_ext")] <- NULL
    params_compare <- params
    params_compare@species_params[c("d", "z_ext")] <- NULL
    expect_unchanged(p2_compare, params_compare)
    
    # supplying ext_mort
    p2 <- setExtMort(params, 3 * params@mu_b)
    comment(p2@mu_b) <- NULL
    expect_equal(p2@mu_b, 3 * params@mu_b)
    
    # has updated time_modified
    expect_false(identical(params@time_modified, p2@time_modified))
})

test_that("setExtMort uses z_ext and d species parameters", {
    params <- NS_params_small
    params@species_params$z0 <- rep(0.2, nrow(params@species_params))
    params@species_params$z_ext <- seq_len(nrow(params@species_params)) / 10
    params@species_params$d <- seq(-0.25, by = 0.01,
                                   length.out = nrow(params@species_params))

    p2 <- setExtMort(params)

    expected <- outer(params@species_params$d, params@w,
                      function(d, w) w^d)
    expected <- sweep(expected, 1, params@species_params$z_ext, "*")
    expected <- sweep(expected, 1, params@species_params$z0, "+")
    expect_equal(p2@mu_b, expected, ignore_attr = TRUE)
})

test_that("setExtMort defaults z_ext to 0 and d to n - 1", {
    params <- NS_params_small
    params@species_params$z_ext <- NULL
    params@species_params$d <- NULL

    p2 <- setExtMort(params, reset = TRUE)

    expect_equal(p2@mu_b,
                 outer(p2@species_params$z0, rep(1, length(p2@w))),
                 ignore_attr = TRUE)
    expect_equal(p2@species_params$z_ext, rep(0, nrow(p2@species_params)))
    expect_equal(p2@species_params$d, p2@species_params$n - 1)
})

test_that("Comment works on mu_b", {
    params <- NS_params_small
    # if no comment, it is set automatically
    ext_mort <- params@mu_b
    params <- setExtMort(params, ext_mort = ext_mort)
    expect_identical(comment(params@mu_b), "set manually")
    
    # comment is stored
    comment(ext_mort) <- "test"
    params <- setExtMort(params, ext_mort = ext_mort)
    expect_identical(comment(params@mu_b), "test")
    
    # if no comment, previous comment is kept
    comment(ext_mort) <- NULL
    params <- setExtMort(params, ext_mort = ext_mort)
    expect_identical(comment(params@mu_b), "test")
    
    # no message when nothing changes
    expect_message(setExtMort(params), NA)
    # but message when a change is not stored due to comment
    params@species_params$z0 <- 1
    expect_message(setExtMort(params),  "has been commented")
    # Can reset
    p <- setExtMort(params, reset = TRUE)
    expect_equal(p@mu_b[1, 1], 1, ignore_attr = TRUE)
    expect_warning(setExtMort(params, ext_mort = ext_mort,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

# second-order bin averaging ----
test_that("setExtMort bin_average path matches the analytic closed form", {
    params <- NS_params_small
    params@species_params$z0 <- rep(0.2, nrow(params@species_params))
    params@species_params$z_ext <- seq_len(nrow(params@species_params)) / 10
    params@species_params$d <- seq(-1, by = 0.3,
                                   length.out = nrow(params@species_params))
    # First species has d == -1 to exercise the logarithmic branch

    params@second_order_w[["bin_average"]] <- TRUE
    p2 <- setExtMort(params)

    # Build the expected bin averages from the closed form
    w <- params@w
    w_next <- params@w + params@dw
    dw <- params@dw
    expected <- matrix(0, nrow = nrow(params@species_params), ncol = length(w))
    for (i in seq_len(nrow(params@species_params))) {
        d <- params@species_params$d[i]
        if (d == -1) {
            avg <- log(w_next / w) / dw
        } else {
            avg <- (w_next^(d + 1) - w^(d + 1)) / ((d + 1) * dw)
        }
        expected[i, ] <- params@species_params$z0[i] +
            params@species_params$z_ext[i] * avg
    }
    expect_equal(p2@mu_b, expected, ignore_attr = TRUE)
})

test_that("setExtMort default (point-sampling) path is unchanged", {
    params <- NS_params_small
    params@species_params$z0 <- rep(0.2, nrow(params@species_params))
    params@species_params$z_ext <- seq_len(nrow(params@species_params)) / 10
    params@species_params$d <- seq(-0.25, by = 0.01,
                                   length.out = nrow(params@species_params))

    # Default flag is FALSE
    expect_false(isTRUE(params@second_order_w[["bin_average"]]))
    p_point <- setExtMort(params)

    expected <- outer(params@species_params$d, params@w,
                      function(d, w) w^d)
    expected <- sweep(expected, 1, params@species_params$z_ext, "*")
    expected <- sweep(expected, 1, params@species_params$z0, "+")
    expect_equal(p_point@mu_b, expected, ignore_attr = TRUE)

    # bin_average differs from point sampling when z_ext != 0
    params@second_order_w[["bin_average"]] <- TRUE
    p_bin <- setExtMort(params)
    expect_true(different(p_bin@mu_b, p_point@mu_b))
})

test_that("setExtMort bin_average leaves constant mortality unchanged", {
    params <- NS_params_small
    params@species_params$z_ext <- 0
    p_point <- setExtMort(params, reset = TRUE)
    params@second_order_w[["bin_average"]] <- TRUE
    p_bin <- setExtMort(params, reset = TRUE)
    expect_equal(p_bin@mu_b, p_point@mu_b, ignore_attr = TRUE)
})

test_that("setExtMort bin_average persists through setParams and projects", {
    params <- NS_params_small
    params@species_params$z_ext <- seq_len(nrow(params@species_params)) / 10
    second_order_w(params) <- c(bin_average = TRUE)
    expect_true(params@second_order_w[["bin_average"]])

    # mu_b is the bin-averaged version after the setParams pipeline
    p_ref <- params
    p_ref@second_order_w[["bin_average"]] <- TRUE
    p_ref <- setExtMort(p_ref, reset = TRUE)
    expect_equal(params@mu_b, p_ref@mu_b, ignore_attr = TRUE)

    sim <- project(params, t_max = 1, t_save = 1)
    expect_true(all(is.finite(sim@n)))
})

test_that("power_law_bin_average matches numerical integration", {
    w <- c(1, 2, 4)
    dw <- c(1, 2, 4)
    for (d in c(-1, -0.5, 0, 0.75, 2)) {
        avg <- power_law_bin_average(w, dw, d)
        num <- vapply(seq_along(w), function(j) {
            integrate(function(x) x^d, w[j], w[j] + dw[j])$value / dw[j]
        }, numeric(1))
        expect_equal(avg, num)
    }
})

# getExtMort ----
test_that("getExtMort works", {
    expect_true(is.ArraySpeciesBySize(getExtMort(NS_params_small)))
    expect_equal(getExtMort(NS_params_small), NS_params_small@mu_b, ignore_attr = TRUE)
})

test_that("Can get and set slot", {
    params <- NS_params_small
    ext_mort <- getExtMort(params)
    expect_identical(ext_mort(params), ext_mort)
    new <- 2 * ext_mort
    comment(new) <- "test"
    ext_mort(params) <- new
    expect_equal(ext_mort(params), new, ignore_attr = TRUE)
    expect_identical(comment(params@mu_b), "test")
})

test_that("setExtMort validates dimensions and deprecated z0 argument", {
    expect_error(setExtMort(NS_params_small, array(0, dim = c(1, 1))))
    expect_warning(p <- setExtMort(NS_params_small, z0 = NS_params_small@mu_b),
                   "deprecated")
    expect_equal(p@mu_b, NS_params_small@mu_b, ignore_attr = TRUE)
})
