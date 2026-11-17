# Tests for sizeIntegral() and its helpers (issue #494)

# helpers ----
test_that("size_dim_labels labels the dimensions of the usual arrays", {
    p <- NS_params_small
    no_sp <- nrow(p@species_params)
    no_w <- length(p@w)
    expect_identical(size_dim_labels(1, "weighting", no_sp, no_w), character(0))
    expect_identical(size_dim_labels(p@w, "weighting", no_sp, no_w), "w")
    expect_identical(size_dim_labels(p@initial_n, "n", no_sp, no_w),
                     c("sp", "w"))
    expect_identical(size_dim_labels(getFMortGear(p), "weighting", no_sp, no_w),
                     c("gear", "sp", "w"))
    # An unnamed extra dimension gets a label of its own
    A <- array(0, dim = c(2, no_sp, no_w))
    expect_identical(size_dim_labels(A, "weighting", no_sp, no_w),
                     c("weighting.dim1", "sp", "w"))
    # "species" and "size" are normalised to mizer's "sp" and "w"
    B <- array(0, dim = c(no_sp, no_w),
               dimnames = list(species = NULL, size = NULL))
    expect_identical(size_dim_labels(B, "weighting", no_sp, no_w), c("sp", "w"))
})

test_that("size_dim_labels rejects arrays of the wrong shape", {
    p <- NS_params_small
    no_sp <- nrow(p@species_params)
    no_w <- length(p@w)
    expect_error(size_dim_labels(1:3, "weighting", no_sp, no_w),
                 "last dimension of `weighting` must run over")
    expect_error(size_dim_labels(array(0, dim = c(2, no_w)), "weighting",
                                 no_sp, no_w),
                 "second-to-last dimension of `weighting` must run over")
    expect_error(size_dim_labels("a", "weighting", no_sp, no_w),
                 "must be numeric")
})

test_that("merge_dim_labels interleaves labels", {
    expect_identical(merge_dim_labels(c("sp", "w"), c("sp", "w")),
                     c("sp", "w"))
    expect_identical(merge_dim_labels(c("sp", "w"), c("gear", "sp", "w")),
                     c("gear", "sp", "w"))
    expect_identical(merge_dim_labels(c("time", "sp", "w"),
                                      c("gear", "sp", "w")),
                     c("time", "gear", "sp", "w"))
    # A shared "time" label is not multiplied out
    expect_identical(merge_dim_labels(c("time", "sp", "w"),
                                      c("time", "gear", "sp", "w")),
                     c("time", "gear", "sp", "w"))
    expect_identical(merge_dim_labels(character(0), c("sp", "w")),
                     c("sp", "w"))
})

test_that("broadcast_dims replicates and permutes correctly", {
    extent <- c(gear = 2L, sp = 3L, w = 4L)
    x <- matrix(1:12, nrow = 3)  # sp x w
    b <- broadcast_dims(x, c("sp", "w"), c("gear", "sp", "w"), extent)
    expect_identical(dim(b), c(2L, 3L, 4L))
    expect_identical(b[1, , ], x)
    expect_identical(b[2, , ], x)
    # a scalar
    s <- broadcast_dims(2, character(0), c("gear", "sp", "w"), extent)
    expect_true(all(s == 2))
    expect_identical(dim(s), c(2L, 3L, 4L))
    # dimensions given in a different order are permuted
    y <- matrix(1:8, nrow = 2)  # gear x w
    b2 <- broadcast_dims(y, c("gear", "w"), c("gear", "sp", "w"), extent)
    expect_identical(b2[, 1, ], y)
    expect_identical(b2[, 3, ], y)
})

test_that("dim_extents catches inconsistent dimensions", {
    expect_error(dim_extents(list(matrix(0, 2, 3), matrix(0, 2, 4)),
                             list(c("sp", "w"), c("sp", "w"))),
                 "'w' dimensions do not have the same length")
    expect_identical(dim_extents(list(matrix(0, 2, 3)), list(c("sp", "w"))),
                     c(sp = 2L, w = 3L))
})

# sizeIntegral ----
test_that("sizeIntegral reproduces the built-in summary functions", {
    p <- NS_params_small
    expect_identical(sizeIntegral(p, weighting = p@w), getBiomass(p))
    expect_identical(sizeIntegral(p), getN(p))
    expect_identical(sizeIntegral(p, weighting = p@w, min_w = 10, max_w = 1000),
                     getBiomass(p, min_w = 10, max_w = 1000))
    expect_equal(sizeIntegral(p, weighting = sweep(p@maturity, 2, p@w, "*")),
                 getSSB(p))
    f <- getFMort(p, drop = FALSE)
    expect_equal(sizeIntegral(p, weighting = sweep(f, 2, p@w, "*")),
                 getYield(p))
    fg <- getFMortGear(p)
    expect_equal(sizeIntegral(p, weighting = sweep(fg, 3, p@w, "*")),
                 getYieldGear(p))
})

test_that("sizeIntegral integrates against a supplied abundance", {
    p <- NS_params_small
    n <- initialN(p)
    expect_identical(sizeIntegral(p, weighting = p@w, n = 2 * n),
                     2 * getBiomass(p))
    # A single number for the weighting factor integrates the abundance density
    expect_identical(sizeIntegral(p, weighting = 1), getN(p))
})

test_that("sizeIntegral gets the shape of the result right", {
    p <- NS_params_small
    sim <- NS_sim_small
    no_sp <- nrow(p@species_params)
    no_gear <- length(dimnames(getFMortGear(p))$gear)

    # MizerParams, weighting over species and size -> vector over species
    b <- sizeIntegral(p, weighting = p@w)
    expect_null(dim(b))
    expect_identical(names(b), as.character(p@species_params$species))

    # MizerSim -> ArrayTimeBySpecies
    bs <- sizeIntegral(sim, weighting = p@w, value_name = "Biomass", units = "g")
    expect_true(is.ArrayTimeBySpecies(bs))
    expect_identical(dim(bs), c(length(getTimes(sim)), no_sp))
    expect_identical(attr(bs, "value_name"), "Biomass")
    expect_identical(attr(bs, "units"), "g")

    # An extra dimension of the weighting is carried through
    yg <- sizeIntegral(p, weighting = getFMortGear(p))
    expect_identical(dim(yg), c(no_gear, no_sp))
    expect_identical(names(dimnames(yg)), c("gear", "sp"))

    # ... also when the abundance has a time dimension
    ygs <- sizeIntegral(sim, weighting = getFMortGear(sim))
    expect_identical(dim(ygs), c(length(getTimes(sim)), no_gear, no_sp))
    expect_identical(names(dimnames(ygs)), c("time", "gear", "sp"))
    # and the time dimension is matched, not multiplied out
    expect_equal(ygs[1, , ], sizeIntegral(p, weighting = getFMortGear(sim)[1, , , ],
                                          n = sim@n[1, , ]))
})

test_that("sizeIntegral respects the size range", {
    p <- NS_params_small
    full <- sizeIntegral(p, weighting = p@w)
    part <- sizeIntegral(p, weighting = p@w, min_w = 10)
    expect_true(all(part < full))
    # min_l is converted to a weight with the species' a and b
    sp <- p@species_params
    expect_equal(sizeIntegral(p, weighting = p@w, min_l = 10),
                 sizeIntegral(p, weighting = p@w, min_w = sp$a * 10^sp$b))
    # An empty size range gives zero
    expect_equal(unname(sizeIntegral(p, weighting = p@w,
                                     min_w = 2 * max(p@w))),
                 rep(0, nrow(sp)))
})

test_that("sizeIntegral rejects invalid input", {
    p <- NS_params_small
    expect_error(sizeIntegral(1), "must be a MizerParams or a MizerSim object")
    expect_error(sizeIntegral(p, weighting = 1:3),
                 "last dimension of `weighting` must run over")
    expect_error(sizeIntegral(p, n = 1:3), "last dimension of `n` must run over")
    expect_error(sizeIntegral(p, n = p@w),
                 "`n` must have a dimension running over the species")
})

# The quadrature scheme ----
test_that("sizeIntegral is unchanged by the default scheme", {
    p <- NS_params_small
    # Explicitly first order: the weighting factor is used as given
    K <- p@w
    expect_identical(unname(sizeIntegral(p, weighting = K)),
                     unname(drop(initialN(p) %*% (K * p@dw))))
})

test_that("sizeIntegral bin-averages the whole weighting factor when asked to", {
    p <- NS_params_small
    second_order_w(p) <- c(bin_average = TRUE)
    # The product maturity * w is averaged as a single weighting factor, which is not
    # the same as averaging the two factors separately.
    K <- sweep(p@maturity, 2, p@w, "*")
    expected <- drop(rowSums(initialN(p) *
                                 sweep(trapezoidal_bin_average(K), 2, p@dw, "*")))
    expect_equal(unname(sizeIntegral(p, weighting = K)), unname(expected))
    wrong <- drop(rowSums(initialN(p) *
        sweep(trapezoidal_bin_average(p@maturity) *
                  rep(trapezoidal_bin_average(p@w), each = nrow(p@maturity)),
              2, p@dw, "*")))
    expect_false(isTRUE(all.equal(unname(expected), unname(wrong))))
})

test_that("sizeIntegral bin-averages the size range mask with the weighting", {
    p <- NS_params_small
    p2 <- p
    second_order_w(p2) <- c(bin_average = TRUE)
    # With a restricted size range the bin straddling the boundary contributes
    # only partially, so the second-order value differs from the first-order
    # one even for the constant weighting of getN().
    expect_false(isTRUE(all.equal(unname(sizeIntegral(p, min_w = 10)),
                                  unname(sizeIntegral(p2, min_w = 10)))))
    # Over the full size range the mask is all ones and nothing changes.
    expect_equal(unname(sizeIntegral(p)), unname(sizeIntegral(p2)))
})

test_that("the two schemes converge as the grid is refined", {
    sp <- NS_species_params_small
    diffs <- vapply(c(50, 100, 200), function(no_w) {
        p1 <- suppressMessages(
            newMultispeciesParams(sp, inter_small, no_w = no_w, info_level = 0))
        p2 <- p1
        second_order_w(p2) <- c(bin_average = TRUE)
        b1 <- sizeIntegral(p1, weighting = p1@w)
        b2 <- sizeIntegral(p2, weighting = p2@w)
        max(abs(b2 - b1) / b1)
    }, numeric(1))
    expect_true(all(diff(diffs) < 0))
})

# The ends of the size range ----
# NS_params_small is of edition 1. Recording another edition in it changes
# nothing else, so these variants differ only in how the size range is bounded.
bounds_params <- function(edition, bin_average = FALSE) {
    p <- NS_params_small
    if (bin_average) second_order_w(p) <- c(bin_average = TRUE)
    defaults_edition(p) <- edition
    p
}
# Limits that lie inside size bins of NS_params_small, well below its top bin
bounds_a <- 3.3
bounds_b <- 1234

test_that("exact bounds are the default from edition 3 on", {
    for (bin_average in c(FALSE, TRUE)) {
        p1 <- bounds_params(1, bin_average)
        p2 <- bounds_params(2, bin_average)
        p3 <- bounds_params(3, bin_average)
        legacy <- sizeIntegral(p1, weighting = p1@w,
                               min_w = bounds_a, max_w = bounds_b)
        exact <- sizeIntegral(p1, weighting = p1@w,
                              min_w = bounds_a, max_w = bounds_b,
                              exact_bounds = TRUE)
        expect_false(isTRUE(all.equal(legacy, exact)))
        expect_identical(sizeIntegral(p2, weighting = p2@w, min_w = bounds_a,
                                      max_w = bounds_b), legacy)
        expect_identical(sizeIntegral(p3, weighting = p3@w, min_w = bounds_a,
                                      max_w = bounds_b), exact)
        expect_identical(sizeIntegral(p3, weighting = p3@w, min_w = bounds_a,
                                      max_w = bounds_b, exact_bounds = FALSE),
                         legacy)
    }
})

test_that("editions 1 and 2 keep selecting whole bins by their left edge", {
    for (bin_average in c(FALSE, TRUE)) {
        for (edition in 1:2) {
            p <- bounds_params(edition, bin_average)
            no_w <- length(p@w)
            K <- get_size_range_array(p, min_w = bounds_a, max_w = bounds_b) *
                rep(p@w, each = nrow(p@species_params))
            K <- sweep(bin_average_weight(K, p), 2, p@dw, "*")
            expected <- drop((initialN(p) * K) %*% rep(1, no_w))
            expect_identical(sizeIntegral(p, weighting = p@w, min_w = bounds_a,
                                          max_w = bounds_b), expected)
        }
    }
})

test_that("exact bounds integrate the boundary bins exactly", {
    no_sp <- nrow(NS_params_small@species_params)
    for (bin_average in c(FALSE, TRUE)) {
        p <- bounds_params(3, bin_average)
        w <- p@w
        dw <- p@dw
        n1 <- matrix(1, no_sp, length(w))
        # With a constant density and weighting the integral is the width of
        # the size range, under both schemes.
        expect_equal(unname(sizeIntegral(p, n = n1, min_w = bounds_a,
                                         max_w = bounds_b)),
                     rep(bounds_b - bounds_a, no_sp))
        # Length limits are converted for each species.
        sp <- p@species_params
        expect_equal(unname(sizeIntegral(p, n = n1, min_l = 10, max_l = 40)),
                     unname(sp$a * (40^sp$b - 10^sp$b)))
        # With the weighting w ...
        biomass <- unname(sizeIntegral(p, weighting = w, n = n1,
                                       min_w = bounds_a, max_w = bounds_b))
        if (bin_average) {
            # ... the trapezoidal rule with w interpolated linearly across
            # the bin is exact.
            expect_equal(biomass, rep((bounds_b^2 - bounds_a^2) / 2, no_sp))
        } else {
            # ... the point weight is multiplied by the fraction of the bin
            # that is inside the size range.
            j <- findInterval(bounds_a, w)
            m <- findInterval(bounds_b, w)
            inner <- (j + 1):(m - 1)
            expected <-
                w[j] * dw[j] * (w[j + 1] - bounds_a) / (w[j + 1] - w[j]) +
                sum(w[inner] * dw[inner]) +
                w[m] * dw[m] * (bounds_b - w[m]) / (w[m + 1] - w[m])
            expect_equal(biomass, rep(expected, no_sp))
        }
    }
})

test_that("exact bounds are additive over adjacent size ranges", {
    for (bin_average in c(FALSE, TRUE)) {
        p <- bounds_params(3, bin_average)
        # A split point inside a bin
        split <- 47.1
        whole <- sizeIntegral(p, weighting = p@w, min_w = bounds_a,
                              max_w = bounds_b)
        parts <- sizeIntegral(p, weighting = p@w, min_w = bounds_a,
                              max_w = split) +
            sizeIntegral(p, weighting = p@w, min_w = split, max_w = bounds_b)
        expect_equal(parts, whole)
    }
})

test_that("exact bounds change continuously with the limits", {
    for (bin_average in c(FALSE, TRUE)) {
        p <- bounds_params(3, bin_average)
        node <- p@w[8]
        below <- sizeIntegral(p, weighting = p@w, min_w = node * (1 - 1e-9))
        at <- sizeIntegral(p, weighting = p@w, min_w = node)
        above <- sizeIntegral(p, weighting = p@w, min_w = node * (1 + 1e-9))
        expect_equal(below, at, tolerance = 1e-8)
        expect_equal(above, at, tolerance = 1e-8)
        expect_true(all(below >= at & at >= above))
    }
})

test_that("limits on grid points need no partial bins", {
    for (bin_average in c(FALSE, TRUE)) {
        p <- bounds_params(3, bin_average)
        j <- 5
        m <- 15
        inside <- seq_along(p@w) >= j & seq_along(p@w) < m
        K <- bin_average_weight(rep(p@w, each = nrow(p@species_params)) *
                                    matrix(1, nrow(p@species_params),
                                           length(p@w)), p)
        expected <- rowSums(initialN(p)[, inside] *
                                sweep(K, 2, p@dw, "*")[, inside])
        expect_equal(sizeIntegral(p, weighting = p@w, min_w = p@w[j],
                                  max_w = p@w[m]),
                     expected)
    }
})

test_that("bins inside the size range get bit-identical weights", {
    for (bin_average in c(FALSE, TRUE)) {
        p <- bounds_params(1, bin_average)
        exact <- size_integral_weights(p, weighting = p@w, min_w = bounds_a,
                                       max_w = bounds_b, exact_bounds = TRUE)
        legacy <- size_integral_weights(p, weighting = p@w, min_w = bounds_a,
                                        max_w = bounds_b, exact_bounds = FALSE)
        expect_identical(exact$labels, c("sp", "w"))
        expect_identical(legacy$labels, c("sp", "w"))
        whole <- p@w >= bounds_a & c(p@w[-1], Inf) <= bounds_b
        expect_identical(exact$weights[, whole], legacy$weights[, whole])
        # Bins outside the range get no weight
        outside <- c(p@w[-1], Inf) <= bounds_a | p@w >= bounds_b
        expect_true(all(exact$weights[, outside] == 0))
    }
})

test_that("an upper limit at the last grid point covers the last bin", {
    for (bin_average in c(FALSE, TRUE)) {
        p <- bounds_params(3, bin_average)
        full <- sizeIntegral(p, weighting = p@w)
        expect_identical(sizeIntegral(p, weighting = p@w, max_w = max(p@w)),
                         full)
        expect_identical(sizeIntegral(p, weighting = p@w, max_w = Inf), full)
        # Over the whole size range nothing changes with the edition
        expect_identical(full, sizeIntegral(p, weighting = p@w,
                                            exact_bounds = FALSE))
    }
})

test_that("exact_bounds is validated", {
    expect_error(sizeIntegral(NS_params_small, exact_bounds = "yes"),
                 "`exact_bounds` must be TRUE, FALSE or NULL")
    expect_error(sizeIntegral(NS_params_small, exact_bounds = NA),
                 "`exact_bounds` must be TRUE, FALSE or NULL")
})

test_that("size_integral_weights builds weights on the resource grid", {
    p <- bounds_params(3)
    lw <- resource_length_params(p)
    # Limits inside bins of the resource grid
    a <- 1.7e-5
    b <- 0.31
    exact <- size_integral_weights(p, min_w = a, max_w = b, resource = TRUE)
    expect_identical(dim(exact$weights), c(1L, length(p@w_full)))
    expect_identical(dimnames(exact$weights)$sp, "Resource")
    expect_equal(sum(exact$weights), b - a)
    # Lengths are converted with the resource's weight-length relationship
    l <- c(0.03, 0.5)
    by_length <- size_integral_weights(p, min_l = l[1], max_l = l[2],
                                       resource = TRUE)
    expect_equal(sum(by_length$weights), lw$a * (l[2]^lw$b - l[1]^lw$b))
    # The left-edge rule counts the bins whose left edge is in the range
    legacy <- size_integral_weights(p, min_w = a, max_w = b, resource = TRUE,
                                    exact_bounds = FALSE)
    in_range <- p@w_full >= a & p@w_full <= b
    expect_equal(sum(legacy$weights), sum(p@dw_full[in_range]))
    # A weighting over the resource grid is applied
    biomass <- size_integral_weights(p, weighting = p@w_full,
                                     min_w = a, max_w = b, resource = TRUE)
    expect_equal(sum(biomass$weights * initialNResource(p)),
                 sum(exact$weights * p@w_full * initialNResource(p)))
})

# The summary functions agree under exact bounds ----
test_that("summary functions with a size range use exact bounds in edition 3", {
    for (bin_average in c(FALSE, TRUE)) {
        p1 <- bounds_params(1, bin_average)
        p3 <- bounds_params(3, bin_average)
        exact_biomass <- sizeIntegral(p1, weighting = p1@w, min_w = bounds_a,
                                      max_w = bounds_b, exact_bounds = TRUE)
        exact_n <- sizeIntegral(p1, min_w = bounds_a, max_w = bounds_b,
                                exact_bounds = TRUE)
        expect_identical(getBiomass(p3, min_w = bounds_a, max_w = bounds_b),
                         exact_biomass)
        expect_identical(getN(p3, min_w = bounds_a, max_w = bounds_b), exact_n)
        # exact_bounds reaches sizeIntegral() through the summary functions
        expect_identical(getBiomass(p1, min_w = bounds_a, max_w = bounds_b,
                                    exact_bounds = TRUE), exact_biomass)
        expect_identical(getN(p3, min_w = bounds_a, max_w = bounds_b,
                              exact_bounds = FALSE),
                         getN(p1, min_w = bounds_a, max_w = bounds_b))
        # ... also with the biomass cutoff
        p1c <- p1
        p1c@species_params$biomass_cutoff <- bounds_a
        p3c <- p3
        p3c@species_params$biomass_cutoff <- bounds_a
        expect_identical(getBiomass(p3c, use_cutoff = TRUE),
                         getBiomass(p3, min_w = bounds_a))
        expect_identical(getBiomass(p1c, use_cutoff = TRUE,
                                    exact_bounds = TRUE),
                         getBiomass(p3, min_w = bounds_a))
        # The indicators are built from the same integrals
        expect_equal(getMeanWeight(p3, min_w = bounds_a, max_w = bounds_b),
                     sum(exact_biomass) / sum(exact_n))
        threshold <- 77
        large <- getProportionOfLargeFish(p3, min_w = bounds_a,
                                          max_w = bounds_b,
                                          threshold_w = threshold)
        small_biomass <- getBiomass(p3, min_w = bounds_a, max_w = threshold)
        expect_equal(large, 1 - sum(small_biomass) / sum(exact_biomass))
        expect_equal(large,
                     sum(getBiomass(p3, min_w = threshold, max_w = bounds_b)) /
                         sum(exact_biomass))
        # Without a size range nothing changes with the edition
        expect_identical(getBiomass(p3), getBiomass(p1))
        expect_identical(getN(p3), getN(p1))
        expect_identical(getSSB(p3), getSSB(p1))
        expect_identical(getYield(p3), getYield(p1))
        expect_identical(getYieldGear(p3), getYieldGear(p1))
    }
})

test_that("summary functions of a simulation use exact bounds in edition 3", {
    sim <- NS_sim_small
    defaults_edition(sim@params) <- 3
    p <- sim@params
    biomass <- getBiomass(sim, min_w = bounds_a, max_w = bounds_b)
    expect_true(is.ArrayTimeBySpecies(biomass))
    for (t in seq_len(dim(sim@n)[[1]])) {
        expect_equal(unclass(biomass)[t, ],
                     sizeIntegral(p, weighting = p@w, n = sim@n[t, , ],
                                  min_w = bounds_a, max_w = bounds_b))
    }
    legacy <- getBiomass(NS_sim_small, min_w = bounds_a, max_w = bounds_b)
    expect_false(isTRUE(all.equal(unclass(biomass), unclass(legacy))))
    # Without a size range nothing changes with the edition
    expect_identical(as.numeric(unclass(getYield(sim))),
                     as.numeric(unclass(getYield(NS_sim_small))))
})
