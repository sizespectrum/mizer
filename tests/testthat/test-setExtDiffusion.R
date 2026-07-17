test_that("setExtDiffusion sets and returns the ext_diffusion array", {
    params <- single_sp_params
    new <- ext_diffusion(params) + 1

    updated <- setExtDiffusion(params, ext_diffusion = new)

    expect_true(is.ArraySpeciesBySize(ext_diffusion(updated)))
    expect_equal(ext_diffusion(updated), new, ignore_attr = TRUE)
    expect_equal(updated@ext_diffusion, new, ignore_attr = TRUE)
})

test_that("setExtDiffusion preserves comments from the supplied array", {
    params <- single_sp_params
    new <- ext_diffusion(params)
    comment(new) <- "custom"

    updated <- setExtDiffusion(params, ext_diffusion = new)

    expect_identical(comment(updated@ext_diffusion), "custom")
})

test_that("setExtDiffusion works", {
    params <- NS_params_small

    # Without ext_diffusion argument, recalculates from D_ext (default 0) and n
    p2 <- setExtDiffusion(params)
    zero_diffusion <- p2@ext_diffusion
    zero_diffusion[] <- 0
    expect_identical(zero_diffusion, p2@ext_diffusion)

    # supplying ext_diffusion
    p2 <- setExtDiffusion(params, 3 * params@ext_diffusion + 1)
    expect_equal(p2@ext_diffusion, 3 * params@ext_diffusion + 1,
                 ignore_attr = TRUE)

    # only ext_diffusion changed
    p2@ext_diffusion <- params@ext_diffusion
    expect_unchanged(p2, params)

    # has updated time_modified
    expect_false(identical(params@time_modified, p2@time_modified))
})

test_that("setExtDiffusion defaults D_ext to 0", {
    # setExtDiffusion owns the D_ext default, so it must supply it even when
    # called standalone, without going through setParams()/validParams().
    params <- NS_params_small
    params@species_params$D_ext <- NULL

    p2 <- setExtDiffusion(params, reset = TRUE)

    expect_equal(p2@species_params$D_ext, rep(0, nrow(p2@species_params)),
                 ignore_attr = TRUE)
    expect_equal(p2@ext_diffusion, p2@ext_diffusion * 0, ignore_attr = TRUE)
})

test_that("setExtDiffusion uses D_ext species param", {
    params <- NS_params_small
    species_params(params)$D_ext <- 0.1
    p2 <- setExtDiffusion(params, reset = TRUE)
    expected <- sweep(outer(species_params(params)[["n"]],
                            w(params), function(x, y) y^(x + 1)),
                      1, species_params(params)[["D_ext"]], "*")
    expect_equal(p2@ext_diffusion, expected, ignore_attr = TRUE)
})

test_that("reset works on ext_diffusion", {
    params <- NS_params_small
    # Set a custom ext_diffusion with a comment
    custom <- params@ext_diffusion
    custom[] <- 1
    comment(custom) <- "custom"
    params <- setExtDiffusion(params, ext_diffusion = custom)
    expect_identical(comment(params@ext_diffusion), "custom")

    # reset = TRUE ignores comment and recalculates from species params (D_ext=0)
    p2 <- setExtDiffusion(params, reset = TRUE)
    expect_null(comment(p2@ext_diffusion))
    expect_equal(p2@ext_diffusion, params@ext_diffusion * 0, ignore_attr = TRUE)
})

test_that("Comment works on ext_diffusion", {
    params <- NS_params_small
    ext_diffusion <- params@ext_diffusion
    # comment is stored
    comment(ext_diffusion) <- "test"
    params <- setExtDiffusion(params, ext_diffusion = ext_diffusion)
    expect_identical(comment(params@ext_diffusion), "test")

    # if no comment, previous comment is kept
    comment(ext_diffusion) <- NULL
    params <- setExtDiffusion(params, ext_diffusion = ext_diffusion)
    expect_identical(comment(params@ext_diffusion), "test")
})

test_that("ext_diffusion() returns ArraySpeciesBySize", {
    expect_true(is.ArraySpeciesBySize(ext_diffusion(NS_params_small)))
    expect_equal(ext_diffusion(NS_params_small), NS_params_small@ext_diffusion,
                 ignore_attr = TRUE)
})

test_that("setExtDiffusion validates dimensions", {
    expect_error(setExtDiffusion(NS_params_small, array(0, dim = c(1, 1))))
})

test_that("Can get and set ext_diffusion slot", {
    params <- NS_params_small
    ed <- ext_diffusion(params)
    new <- ed + 1
    comment(new) <- "test"
    ext_diffusion(params) <- new
    expect_equal(ext_diffusion(params), new, ignore_attr = TRUE)
    expect_identical(comment(params@ext_diffusion), "test")
})

test_that("setExtDiffusion bin_average path matches the analytic bin average", {
    params <- NS_params_small
    params@species_params$D_ext <- seq_len(nrow(params@species_params)) / 10
    params@species_params$n <- seq(0.6, by = 0.1,
                                   length.out = nrow(params@species_params))

    params@second_order_w[["bin_average"]] <- TRUE
    p2 <- setExtDiffusion(params)

    # Expected: D_ext times the exact bin average of w^(n+1) over each bin
    w <- params@w
    w_next <- params@w + params@dw
    dw <- params@dw
    expected <- matrix(0, nrow = nrow(params@species_params), ncol = length(w))
    for (i in seq_len(nrow(params@species_params))) {
        e <- params@species_params$n[i] + 1            # exponent of w
        avg <- (w_next^(e + 1) - w^(e + 1)) / ((e + 1) * dw)
        expected[i, ] <- params@species_params$D_ext[i] * avg
    }
    expect_equal(p2@ext_diffusion, expected, ignore_attr = TRUE)

    # Differs from the point-sampled default
    params@second_order_w[["bin_average"]] <- FALSE
    p1 <- setExtDiffusion(params)
    expect_false(isTRUE(all.equal(p1@ext_diffusion, p2@ext_diffusion)))
})
