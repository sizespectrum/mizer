test_that("setExtDiffusion works", {
    params <- NS_params

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

test_that("setExtDiffusion uses D_ext species param", {
    params <- NS_params
    species_params(params)$D_ext <- 0.1
    p2 <- setExtDiffusion(params, reset = TRUE)
    expected <- sweep(outer(species_params(params)[["n"]],
                            w(params), function(x, y) y^(x + 1)),
                      1, species_params(params)[["D_ext"]], "*")
    expect_equal(p2@ext_diffusion, expected, ignore_attr = TRUE)
})

test_that("reset works on ext_diffusion", {
    params <- NS_params
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
    params <- NS_params
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
    expect_true(is.ArraySpeciesBySize(ext_diffusion(NS_params)))
    expect_equal(ext_diffusion(NS_params), NS_params@ext_diffusion,
                 ignore_attr = TRUE)
})

test_that("setExtDiffusion validates dimensions", {
    expect_error(setExtDiffusion(NS_params, array(0, dim = c(1, 1))))
})

test_that("Can get and set ext_diffusion slot", {
    params <- NS_params
    ed <- ext_diffusion(params)
    new <- ed + 1
    comment(new) <- "test"
    ext_diffusion(params) <- new
    expect_equal(ext_diffusion(params), new, ignore_attr = TRUE)
    expect_identical(comment(params@ext_diffusion), "test")
})
