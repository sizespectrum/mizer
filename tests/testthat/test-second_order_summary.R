# Tests for the second-order (bin-averaged) summary integrals (issue #380)

# trapezoidal_bin_average helper ----
test_that("trapezoidal_bin_average trapezoid-averages a vector", {
    K <- c(2, 4, 8, 16)
    expect_equal(trapezoidal_bin_average(K), c(3, 6, 12, 16))
})

test_that("trapezoidal_bin_average averages along the last dimension of an array", {
    # matrix (rows = species, cols = size)
    M <- rbind(c(1, 3, 5), c(2, 4, 8))
    expect_equal(trapezoidal_bin_average(M),
                 rbind(c(2, 4, 5), c(3, 6, 8)))
    # 3D array, average over the last (size) axis
    A <- array(1:12, dim = c(2, 2, 3))
    Abar <- trapezoidal_bin_average(A)
    expect_equal(dim(Abar), dim(A))
    # last bin unchanged (one-sided)
    expect_equal(Abar[, , 3], A[, , 3])
    # interior bins are the trapezoid average
    expect_equal(Abar[, , 1], 0.5 * (A[, , 1] + A[, , 2]))
})

test_that("trapezoidal_bin_average is exact for the monomial weight w", {
    p <- NS_params_small
    w <- p@w
    dw <- p@dw
    no_w <- length(w)
    wbar <- trapezoidal_bin_average(w)
    # interior bins equal (w_{j+1}^2 - w_j^2) / (2 dw_j)
    w_edge <- c(w, w[no_w] + dw[no_w])
    exact <- (w_edge[-1]^2 - w_edge[-(no_w + 1)]^2) / (2 * dw)
    expect_equal(wbar[-no_w], exact[-no_w])
})

test_that("bin_average_weight is gated on second_order_w", {
    p <- NS_params_small
    K <- p@w
    # default: unchanged
    expect_identical(bin_average_weight(K, p), K)
    # second-order: bin-averaged
    second_order_w(p) <- c(bin_average = TRUE)
    expect_identical(bin_average_weight(K, p), trapezoidal_bin_average(K))
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
    weight <- sweep(trapezoidal_bin_average(K), 2, p@dw, "*")
    expect_equal(unname(getSSB(p)),
                 unname(rowSums(p@initial_n * weight)))
})

test_that("second-order getMeanLength uses the bin-averaged length weight", {
    p <- NS_params_small
    second_order_w(p) <- c(bin_average = TRUE)
    l <- length_at_size(p)
    weight <- sweep(trapezoidal_bin_average(l), 2, p@dw, "*")
    expect_equal(getMeanLength(p),
                 sum(p@initial_n * weight) / sum(getN(p)))
    # and it moves the answer away from the default scheme
    expect_false(isTRUE(all.equal(getMeanLength(p),
                                  getMeanLength(NS_params_small))))
})

test_that("second-order getYield uses the bin-averaged F*w weight", {
    p <- NS_params_small
    second_order_w(p) <- c(bin_average = TRUE)
    f <- getFMort(p, drop = FALSE)
    K <- sweep(f, 2, p@w, "*")
    weight <- sweep(trapezoidal_bin_average(K), 2, p@dw, "*")
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


# Diet and trophic level must use the encounter quadrature (issue #474) ----

# When second-order bin-averaging is on, the prey-bin integral is folded into
# the Fourier-transformed kernel by setPredKernel(). A summary function that
# also bin-averages its prey weight applies that quadrature twice and inflates
# the result by (1 + beta) / 2. These tests pin the identities that catch it.

custom_kernel_params <- function(params) {
    pk <- pred_kernel(params)
    comment(pk) <- "set manually"
    setPredKernel(params, pred_kernel = pk)
}

test_that("encounter_kernel reproduces the encounter rate in both modes", {
    check <- function(params) {
        no_w <- length(params@w)
        no_w_full <- length(params@w_full)
        idx_sp <- (no_w_full - no_w + 1):no_w_full
        n <- initialN(params)
        n_pp <- initialNResource(params)
        kernel <- encounter_kernel(params)
        n_eff <- sweep(params@interaction %*% n, 2, params@w * params@dw, "*")
        species <- rowSums(sweep(kernel[, , idx_sp, drop = FALSE], c(1, 3),
                                 n_eff, "*"), dims = 2)
        resource <- params@species_params$interaction_resource *
            rowSums(sweep(kernel, 3, params@w_full * params@dw_full * n_pp,
                          "*"), dims = 2)
        direct <- params@search_vol * (species + resource) +
            params@ext_encounter
        expect_equal(as.vector(direct), as.vector(getEncounter(params)),
                     tolerance = 1e-4)
    }
    p <- NS_params_small
    check(p)
    check(custom_kernel_params(p))
    second_order_w(p) <- c(bin_average = TRUE)
    check(p)
    check(custom_kernel_params(p))
})

test_that("getDiet summed over prey equals the consumption rate in both modes", {
    check <- function(params) {
        total <- rowSums(getDiet(params, proportion = FALSE), dims = 2)
        consumption <- getEncounter(params) * (1 - getFeedingLevel(params))
        # Outside a species' size range the diet is set to zero, so compare
        # only where the abundance is positive.
        mask <- initialN(params) > 0
        expect_equal(as.vector(total)[mask], as.vector(consumption)[mask])
    }
    p <- NS_params_small
    check(p)
    check(custom_kernel_params(p))
    second_order_w(p) <- c(bin_average = TRUE)
    check(p)
    check(custom_kernel_params(p))
})

test_that("getTrophicLevel gives 2 for a predator whose prey all have level 1", {
    # At the smallest size on the grid no consumer has yet been assigned a
    # trophic level above 1, and choosing w_R above the whole grid makes the
    # resource trophic level 1 as well. The trophic-level-weighted encounter in
    # the numerator then equals the plain encounter in the denominator, so the
    # trophic level must come out as exactly 2. It does so only if numerator and
    # denominator use the same quadrature.
    check <- function(params) {
        w_R <- max(params@w_full) * 10
        tl <- getTrophicLevel(params, w_R = w_R)
        first <- which(params@w_min_idx == 1)
        expect_equal(as.vector(tl[first, 1]), rep(2, length(first)),
                     tolerance = 1e-4)
    }
    p <- NS_params_small
    check(p)
    check(custom_kernel_params(p))
    second_order_w(p) <- c(bin_average = TRUE)
    check(p)
    check(custom_kernel_params(p))
})

# Calibration and matching ----

test_that("the calibration and matching functions integrate in both modes", {
    check <- function(params) {
        cutoff <- params@w[10]
        species_params(params)$biomass_cutoff <- cutoff
        species_params(params)$number_cutoff <- cutoff
        biomass <- getBiomass(params, use_cutoff = TRUE)
        number <- getN(params, min_w = cutoff)
        species_params(params)$biomass_observed <- as.numeric(biomass) *
            c(0.5, 2, 4)
        species_params(params)$number_observed <- as.numeric(number) *
            c(4, 0.5, 2)
        observed_biomass <- species_params(params)$biomass_observed
        observed_number <- species_params(params)$number_observed

        # Calibrating brings the total onto the total of the observations
        expect_equal(sum(getBiomass(calibrateBiomass(params),
                                    use_cutoff = TRUE)),
                     sum(observed_biomass))
        expect_equal(sum(getN(calibrateNumber(params), min_w = cutoff)),
                     sum(observed_number))

        # Matching brings each species onto its own observation
        matched <- suppressMessages(matchBiomasses(params))
        expect_equal(unname(getBiomass(matched, use_cutoff = TRUE)),
                     unname(observed_biomass))
        matched <- suppressMessages(matchNumbers(params))
        expect_equal(unname(getN(matched, min_w = cutoff)),
                     unname(observed_number))
    }
    p <- NS_params_small
    check(p)
    second_order_w(p) <- c(bin_average = TRUE)
    check(p)
})
