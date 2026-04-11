test_that("setDiffusion sets and returns the diffusion array", {
    params <- newSingleSpeciesParams()
    new <- diffusion(params) + 1

    updated <- setDiffusion(params, diffusion = new)

    expect_true(is.ArraySpeciesBySize(diffusion(updated)))
    expect_equal(diffusion(updated), new, ignore_attr = TRUE)
    expect_equal(updated@ext_diffusion, new, ignore_attr = TRUE)
})

test_that("setDiffusion preserves comments from the supplied array", {
    params <- newSingleSpeciesParams()
    new <- diffusion(params)
    comment(new) <- "custom"

    updated <- setDiffusion(params, diffusion = new)

    expect_identical(comment(updated@ext_diffusion), "custom")
})

test_that("setDiffusion validates dimensions and accessors delegate correctly", {
    params <- newSingleSpeciesParams()

    expect_error(setDiffusion(params, diffusion = array(0, dim = c(1, 1))))

    new <- diffusion(params) + 2
    diffusion(params) <- new
    expect_identical(params@ext_diffusion, new)
})
