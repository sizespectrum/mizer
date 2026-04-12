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

# getDiffusion ----
test_that("getDiffusion returns ArraySpeciesBySize with correct dimensions", {
    params <- newSingleSpeciesParams()
    d <- getDiffusion(params)
    expect_true(is.ArraySpeciesBySize(d))
    expect_identical(dim(d), dim(params@initial_n))
})

test_that("getDiffusion includes ext_diffusion", {
    params <- newSingleSpeciesParams()
    d_base <- getDiffusion(params)
    # Adding a constant to ext_diffusion should shift getDiffusion by the same amount
    params@ext_diffusion[] <- 1
    d_with_ext <- getDiffusion(params)
    expect_equal(d_with_ext, d_base + 1, ignore_attr = TRUE)
})

test_that("getDiffusion dispatches via rates_funcs", {
    params <- newSingleSpeciesParams()
    e <- globalenv()
    e$constant_diffusion <- function(params, n, n_pp, n_other, t, feeding_level, ...) {
        array(42, dim = dim(params@initial_n), dimnames = dimnames(params@initial_n))
    }
    params@rates_funcs$Diffusion <- "constant_diffusion"
    d <- getDiffusion(params)
    expect_true(all(d == 42))
})

test_that("r$diffusion is included in getRates output", {
    params <- newSingleSpeciesParams()
    r <- getRates(params)
    expect_true("diffusion" %in% names(r))
    expect_identical(dim(r$diffusion), dim(params@initial_n))
})

test_that("r$diffusion matches getDiffusion", {
    params <- newSingleSpeciesParams()
    r <- getRates(params)
    expect_equal(r$diffusion, getDiffusion(params), ignore_attr = TRUE)
})
