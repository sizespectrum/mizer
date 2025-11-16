test_that("pred kernel functions", {
    ppmr <- NS_params@w_full / NS_params@w_full[1]
    phi <- lognormal_pred_kernel(ppmr, beta = 100, sigma = 2)
    expect_true(all(phi > 0))

    phit <- truncated_lognormal_pred_kernel(ppmr, beta = 100, sigma = 2)
    expect_equal(phi[2], phit[2])
    expect_equal(phit[61], 0)

    phib <- box_pred_kernel(ppmr = 1:5, ppmr_min = 2, ppmr_max = 4)
    expect_identical(phib[1], 0)
    expect_identical(phib[2], 1)
    expect_identical(phib[5], 0)

    phip <- power_law_pred_kernel(1:100, kernel_exp = 0,
                                  kernel_l_l = log(25), kernel_u_l = 1000,
                                  kernel_l_r = log(75), kernel_u_r = 1000)
    expect_equal(phip[10], 0)
    expect_equal(phip[30], 1)
    expect_equal(phip[90], 0)
})

test_that("box_pred_kernel returns correct values", {
    # Test basic box kernel
    ppmr <- 1:10
    result <- box_pred_kernel(ppmr, ppmr_min = 3, ppmr_max = 7)

    # Values below min should be 0
    expect_equal(result[1:2], c(0, 0))
    # Values within range should be 1
    expect_equal(result[3:7], rep(1, 5))
    # Values above max should be 0
    expect_equal(result[8:10], c(0, 0, 0))
})

test_that("box_pred_kernel handles boundary values correctly", {
    ppmr <- c(1.9, 2.0, 2.1, 3.9, 4.0, 4.1)
    result <- box_pred_kernel(ppmr, ppmr_min = 2, ppmr_max = 4)

    # Value just below min should be 0
    expect_equal(result[1], 0)
    # Value at min should be 1
    expect_equal(result[2], 1)
    # Value just above min should be 1
    expect_equal(result[3], 1)
    # Value just below max should be 1
    expect_equal(result[4], 1)
    # Value at max should be 1
    expect_equal(result[5], 1)
    # Value just above max should be 0
    expect_equal(result[6], 0)
})

test_that("box_pred_kernel validates parameters", {
    ppmr <- 1:10
    # ppmr_min must be less than ppmr_max
    expect_error(box_pred_kernel(ppmr, ppmr_min = 5, ppmr_max = 3))
    expect_error(box_pred_kernel(ppmr, ppmr_min = 5, ppmr_max = 5))
})

test_that("box_pred_kernel returns vector of correct length", {
    ppmr <- 1:100
    result <- box_pred_kernel(ppmr, ppmr_min = 10, ppmr_max = 50)
    expect_equal(length(result), length(ppmr))
})

test_that("lognormal_pred_kernel returns correct values", {
    ppmr <- c(0.5, 1, 10, 100, 1000)
    result <- lognormal_pred_kernel(ppmr, beta = 100, sigma = 2)

    # All values should be non-negative
    expect_true(all(result >= 0))
    # Length should match input
    expect_equal(length(result), length(ppmr))
    # Values below 1 should be 0 (predator smaller than prey)
    expect_equal(result[ppmr < 1], rep(0, sum(ppmr < 1)))
    # Peak should be near beta
    ppmr_fine <- seq(1, 200, by = 1)
    result_fine <- lognormal_pred_kernel(ppmr_fine, beta = 100, sigma = 1)
    peak_idx <- which.max(result_fine)
    expect_true(abs(ppmr_fine[peak_idx] - 100) < 5)  # Peak within 5 of beta
})

test_that("truncated_lognormal_pred_kernel truncates correctly", {
    ppmr <- seq(0.5, 41000, by = 100)
    result_ln <- lognormal_pred_kernel(ppmr, beta = 100, sigma = 2)
    result_tln <- truncated_lognormal_pred_kernel(ppmr, beta = 100, sigma = 2)

    # Length should match input
    expect_equal(length(result_tln), length(ppmr))
    # Should be zero at ppmr < 1
    expect_equal(result_tln[ppmr < 1], rep(0, sum(ppmr < 1)))
    # Should be zero at ppmr > exp(log(100) + 3 * 2) = 40342.88
    expect_equal(result_tln[ppmr > 40342.88], rep(0, sum(ppmr > 40342.88)))
    # Should match lognormal in the middle range
    expect_equal(result_ln[50], result_tln[50])
})

test_that("power_law_pred_kernel returns correct shape", {
    ppmr <- 1:100
    result <- power_law_pred_kernel(ppmr, kernel_exp = 0,
                                    kernel_l_l = log(25), kernel_u_l = 1000,
                                    kernel_l_r = log(75), kernel_u_r = 1000)

    # All values should be non-negative
    expect_true(all(result >= 0))
    # Length should match input
    expect_equal(length(result), length(ppmr))
    # Should be close to 0 outside the range
    expect_true(result[10] < 0.1)
    expect_true(result[90] < 0.1)
    # Should be close to 1 in the middle of the range
    expect_true(result[50] > 0.9)
})
