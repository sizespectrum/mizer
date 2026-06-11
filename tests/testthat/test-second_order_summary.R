# Tests for the second-order (bin-averaged) summary integrals (issue #380)

# bin_average_weight helper ----
test_that("bin_average_weight trapezoid-averages a vector", {
    K <- c(2, 4, 8, 16)
    expect_equal(bin_average_weight(K), c(3, 6, 12, 16))
})

test_that("bin_average_weight averages along the last dimension of an array", {
    # matrix (rows = species, cols = size)
    M <- rbind(c(1, 3, 5), c(2, 4, 8))
    expect_equal(bin_average_weight(M),
                 rbind(c(2, 4, 5), c(3, 6, 8)))
    # 3D array, average over the last (size) axis
    A <- array(1:12, dim = c(2, 2, 3))
    Abar <- bin_average_weight(A)
    expect_equal(dim(Abar), dim(A))
    # last bin unchanged (one-sided)
    expect_equal(Abar[, , 3], A[, , 3])
    # interior bins are the trapezoid average
    expect_equal(Abar[, , 1], 0.5 * (A[, , 1] + A[, , 2]))
})

test_that("bin_average_weight is exact for the monomial weight w", {
    p <- NS_params_small
    w <- p@w
    dw <- p@dw
    no_w <- length(w)
    wbar <- bin_average_weight(w)
    # interior bins equal (w_{j+1}^2 - w_j^2) / (2 dw_j)
    w_edge <- c(w, w[no_w] + dw[no_w])
    exact <- (w_edge[-1]^2 - w_edge[-(no_w + 1)]^2) / (2 * dw)
    expect_equal(wbar[-no_w], exact[-no_w])
})

test_that("bin_average_summary_weight is gated on second_order_w", {
    p <- NS_params_small
    K <- p@w
    # default: unchanged
    expect_identical(bin_average_summary_weight(K, p), K)
    # second-order: bin-averaged
    second_order_w(p) <- c(bin_average = TRUE)
    expect_identical(bin_average_summary_weight(K, p), bin_average_weight(K))
})

# Default path is byte-identical ----
test_that("default summary integrals are unchanged (byte-identical)", {
    p <- NS_params_small
    w_dw <- p@w * p@dw
    expect_equal(unname(getBiomass(p)),
                 as.vector(p@initial_n %*% w_dw))
    expect_equal(unname(getSSB(p)),
                 as.vector((p@initial_n * p@maturity) %*% w_dw))
    f <- getFMort(p, drop = FALSE)
    expect_equal(unname(getYield(p)),
                 unname(rowSums(sweep(p@initial_n, 2, w_dw, "*") * f)))
    fg <- getFMortGear(p)
    yg <- apply(sweep(fg, c(2, 3), sweep(p@initial_n, 2, w_dw, "*"), "*"),
                c(1, 2), sum)
    expect_equal(getYieldGear(p), yg)
})

test_that("getN is identical in both modes (zeroth moment already exact)", {
    p <- NS_params_small
    p2 <- p
    second_order_w(p2) <- c(bin_average = TRUE)
    expect_identical(getN(p), getN(p2))
    sim2 <- NS_sim_small
    sim2@params <- p2
    # Compare the numeric values (the attached params metadata differs)
    strip <- function(x) array(as.numeric(x), dim = dim(x))
    expect_identical(strip(getN(NS_sim_small)), strip(getN(sim2)))
})

# Second-order path matches the analytic bin-averaged weight ----
test_that("second-order getBiomass matches the exact (w_{j+1}^2-w_j^2)/2 weight", {
    p <- NS_params_small
    second_order_w(p) <- c(bin_average = TRUE)
    w <- p@w
    dw <- p@dw
    no_w <- length(w)
    wbar <- w
    wbar[-no_w] <- 0.5 * (w[-no_w] + w[-1])
    weight <- wbar * dw
    expect_equal(unname(getBiomass(p)),
                 as.vector(p@initial_n %*% weight))
})

test_that("second-order getSSB uses the bin-averaged maturity*w weight", {
    p <- NS_params_small
    second_order_w(p) <- c(bin_average = TRUE)
    K <- sweep(p@maturity, 2, p@w, "*")
    weight <- sweep(bin_average_weight(K), 2, p@dw, "*")
    expect_equal(unname(getSSB(p)),
                 unname(rowSums(p@initial_n * weight)))
})

test_that("second-order getYield uses the bin-averaged F*w weight", {
    p <- NS_params_small
    second_order_w(p) <- c(bin_average = TRUE)
    f <- getFMort(p, drop = FALSE)
    K <- sweep(f, 2, p@w, "*")
    weight <- sweep(bin_average_weight(K), 2, p@dw, "*")
    expect_equal(unname(getYield(p)),
                 unname(rowSums(weight * p@initial_n)))
})

# Second order moves results in the expected direction on a coarse grid ----
test_that("second-order biomass differs from default on a coarse grid", {
    p <- NS_params_small
    p2 <- p
    second_order_w(p2) <- c(bin_average = TRUE)
    expect_false(isTRUE(all.equal(unname(getBiomass(p)),
                                  unname(getBiomass(p2)))))
})

# Second-order converges to the default under grid refinement ----
test_that("second-order biomass converges to default as the grid is refined", {
    # On a finer grid the (beta+1)/2 left-edge bias shrinks, so the two
    # results agree to higher relative accuracy.
    coarse <- NS_params_small  # 20 bins
    fine <- newMultispeciesParams(NS_species_params_small, inter_small,
                                  no_w = 200)
    rel_diff <- function(params) {
        params2 <- params
        second_order_w(params2) <- c(bin_average = TRUE)
        b1 <- getBiomass(params)
        b2 <- getBiomass(params2)
        max(abs(b2 - b1) / b1)
    }
    expect_lt(rel_diff(fine), rel_diff(coarse))
})
