# Tests for sizeIntegral() and its helpers (issue #494)

# helpers ----
test_that("size_dim_labels labels the dimensions of the usual arrays", {
    p <- NS_params_small
    no_sp <- nrow(p@species_params)
    no_w <- length(p@w)
    expect_identical(size_dim_labels(1, "weight", no_sp, no_w), character(0))
    expect_identical(size_dim_labels(p@w, "weight", no_sp, no_w), "w")
    expect_identical(size_dim_labels(p@initial_n, "n", no_sp, no_w),
                     c("sp", "w"))
    expect_identical(size_dim_labels(getFMortGear(p), "weight", no_sp, no_w),
                     c("gear", "sp", "w"))
    # An unnamed extra dimension gets a label of its own
    A <- array(0, dim = c(2, no_sp, no_w))
    expect_identical(size_dim_labels(A, "weight", no_sp, no_w),
                     c("weight.dim1", "sp", "w"))
    # "species" and "size" are normalised to mizer's "sp" and "w"
    B <- array(0, dim = c(no_sp, no_w),
               dimnames = list(species = NULL, size = NULL))
    expect_identical(size_dim_labels(B, "weight", no_sp, no_w), c("sp", "w"))
})

test_that("size_dim_labels rejects arrays of the wrong shape", {
    p <- NS_params_small
    no_sp <- nrow(p@species_params)
    no_w <- length(p@w)
    expect_error(size_dim_labels(1:3, "weight", no_sp, no_w),
                 "last dimension of `weight` must run over")
    expect_error(size_dim_labels(array(0, dim = c(2, no_w)), "weight",
                                 no_sp, no_w),
                 "second-to-last dimension of `weight` must run over")
    expect_error(size_dim_labels("a", "weight", no_sp, no_w),
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
    expect_identical(sizeIntegral(p, weight = p@w), getBiomass(p))
    expect_identical(sizeIntegral(p), getN(p))
    expect_identical(sizeIntegral(p, weight = p@w, min_w = 10, max_w = 1000),
                     getBiomass(p, min_w = 10, max_w = 1000))
    expect_equal(sizeIntegral(p, weight = sweep(p@maturity, 2, p@w, "*")),
                 getSSB(p))
    f <- getFMort(p, drop = FALSE)
    expect_equal(sizeIntegral(p, weight = sweep(f, 2, p@w, "*")),
                 getYield(p))
    fg <- getFMortGear(p)
    expect_equal(sizeIntegral(p, weight = sweep(fg, 3, p@w, "*")),
                 getYieldGear(p))
})

test_that("sizeIntegral integrates against a supplied abundance", {
    p <- NS_params_small
    n <- initialN(p)
    expect_identical(sizeIntegral(p, weight = p@w, n = 2 * n),
                     2 * getBiomass(p))
    # A single number for the weight integrates the abundance density
    expect_identical(sizeIntegral(p, weight = 1), getN(p))
})

test_that("sizeIntegral gets the shape of the result right", {
    p <- NS_params_small
    sim <- NS_sim_small
    no_sp <- nrow(p@species_params)
    no_gear <- length(dimnames(getFMortGear(p))$gear)

    # MizerParams, weight over species and size -> vector over species
    b <- sizeIntegral(p, weight = p@w)
    expect_null(dim(b))
    expect_identical(names(b), as.character(p@species_params$species))

    # MizerSim -> ArrayTimeBySpecies
    bs <- sizeIntegral(sim, weight = p@w, value_name = "Biomass", units = "g")
    expect_true(is.ArrayTimeBySpecies(bs))
    expect_identical(dim(bs), c(length(getTimes(sim)), no_sp))
    expect_identical(attr(bs, "value_name"), "Biomass")
    expect_identical(attr(bs, "units"), "g")

    # An extra dimension of the weight is carried through
    yg <- sizeIntegral(p, weight = getFMortGear(p))
    expect_identical(dim(yg), c(no_gear, no_sp))
    expect_identical(names(dimnames(yg)), c("gear", "sp"))

    # ... also when the abundance has a time dimension
    ygs <- sizeIntegral(sim, weight = getFMortGear(sim))
    expect_identical(dim(ygs), c(length(getTimes(sim)), no_gear, no_sp))
    expect_identical(names(dimnames(ygs)), c("time", "gear", "sp"))
    # and the time dimension is matched, not multiplied out
    expect_equal(ygs[1, , ], sizeIntegral(p, weight = getFMortGear(sim)[1, , , ],
                                          n = sim@n[1, , ]))
})

test_that("sizeIntegral respects the size range", {
    p <- NS_params_small
    full <- sizeIntegral(p, weight = p@w)
    part <- sizeIntegral(p, weight = p@w, min_w = 10)
    expect_true(all(part < full))
    # min_l is converted to a weight with the species' a and b
    sp <- p@species_params
    expect_equal(sizeIntegral(p, weight = p@w, min_l = 10),
                 sizeIntegral(p, weight = p@w, min_w = sp$a * 10^sp$b))
    # An empty size range gives zero
    expect_equal(unname(sizeIntegral(p, weight = p@w,
                                     min_w = 2 * max(p@w))),
                 rep(0, nrow(sp)))
})

test_that("sizeIntegral rejects invalid input", {
    p <- NS_params_small
    expect_error(sizeIntegral(1), "must be a MizerParams or a MizerSim object")
    expect_error(sizeIntegral(p, weight = 1:3),
                 "last dimension of `weight` must run over")
    expect_error(sizeIntegral(p, n = 1:3), "last dimension of `n` must run over")
    expect_error(sizeIntegral(p, n = p@w),
                 "`n` must have a dimension running over the species")
})

# The quadrature scheme ----
test_that("sizeIntegral is unchanged by the default scheme", {
    p <- NS_params_small
    # Explicitly first order: the weight is used as given
    K <- p@w
    expect_identical(unname(sizeIntegral(p, weight = K)),
                     unname(drop(initialN(p) %*% (K * p@dw))))
})

test_that("sizeIntegral bin-averages the whole weight when asked to", {
    p <- NS_params_small
    second_order_w(p) <- c(bin_average = TRUE)
    # The product maturity * w is averaged as a single weight, which is not
    # the same as averaging the two factors separately.
    K <- sweep(p@maturity, 2, p@w, "*")
    expected <- drop(rowSums(initialN(p) *
                                 sweep(trapezoidal_bin_average(K), 2, p@dw, "*")))
    expect_equal(unname(sizeIntegral(p, weight = K)), unname(expected))
    wrong <- drop(rowSums(initialN(p) *
        sweep(trapezoidal_bin_average(p@maturity) *
                  rep(trapezoidal_bin_average(p@w), each = nrow(p@maturity)),
              2, p@dw, "*")))
    expect_false(isTRUE(all.equal(unname(expected), unname(wrong))))
})

test_that("sizeIntegral bin-averages the size range mask with the weight", {
    p <- NS_params_small
    p2 <- p
    second_order_w(p2) <- c(bin_average = TRUE)
    # With a restricted size range the bin straddling the boundary contributes
    # only partially, so the second-order value differs from the first-order
    # one even for the constant weight of getN().
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
        b1 <- sizeIntegral(p1, weight = p1@w)
        b2 <- sizeIntegral(p2, weight = p2@w)
        max(abs(b2 - b1) / b1)
    }, numeric(1))
    expect_true(all(diff(diffs) < 0))
})
