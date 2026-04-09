test_that("setDiffusion sets and returns the diffusion array", {
    params <- NS_params
    new <- params@diffusion + 1

    updated <- setDiffusion(params, diffusion = new)

    expect_identical(diffusion(updated), new)
    expect_identical(updated@diffusion, new)
})

test_that("setDiffusion preserves comments from the supplied array", {
    params <- NS_params
    new <- params@diffusion
    comment(new) <- "custom"

    updated <- setDiffusion(params, diffusion = new)

    expect_identical(comment(updated@diffusion), "custom")
})

test_that("setDiffusion validates dimensions and accessors delegate correctly", {
    params <- NS_params

    expect_error(setDiffusion(params, diffusion = array(0, dim = c(1, 1))))

    new <- params@diffusion + 2
    diffusion(params) <- new
    expect_identical(params@diffusion, new)
})
