# A single-species model in which the species' rates do not depend on its own
# abundance (a fixed community without cannibalism), under the given flux
# scheme, with or without size diffusion, settled onto its steady state for
# fixed rates. For such a model the steady state for fixed rates is a steady
# state of the model, so R0 = RDI / required recruitment there.
r0_params <- function(flux = "upwind", bin_average = FALSE,
                      diffusion = FALSE) {
    params <- single_sp_params
    second_order_w(params) <- flux
    second_order_w(params) <- c(bin_average = bin_average)
    if (diffusion) {
        n <- params@species_params$n[1]
        ext_diffusion(params)[1, ] <- 0.1 * params@w^(n + 1)
    }
    steadySingleSpecies(params, info_level = 0)
}

# getR0 ----
test_that("getR0() is RDI over required recruitment at a steady state", {
    schemes <- list(list("upwind", FALSE, FALSE),
                    list("upwind", FALSE, TRUE),
                    list("van_leer", TRUE, TRUE),
                    list("centred", TRUE, TRUE))
    for (scheme in schemes) {
        params <- do.call(r0_params, scheme)
        expect_equal(getR0(params), getRDI(params) / getRequiredRDD(params),
                     tolerance = 1e-8)
    }
})

test_that("getR0() calibrates erepro to a steady state of project()", {
    # The case the hand-written upwind recursion gets wrong
    params <- r0_params("centred", TRUE, TRUE)
    params <- setRateFunction(params, "RDD", "noRDD")
    params@species_params$erepro <-
        params@species_params$erepro / getR0(params)
    expect_equal(getR0(params), c("Target species" = 1))

    sim <- project(params, t_max = 1, t_save = 1, progress_bar = FALSE)
    expect_equal(finalN(sim), initialN(params), tolerance = 1e-10,
                 ignore_attr = TRUE)
})

test_that("getR0() is 1 / (1 - r) at a Beverton-Holt steady state", {
    params <- r0_params("van_leer", TRUE, TRUE)
    params <- setBevertonHolt(params, reproduction_level = 0.4)
    expect_equal(getR0(params), c("Target species" = 1 / 0.6), tolerance = 1e-8)

    params <- NS_params_steady_small
    expect_equal(getR0(params), 1 / (1 - reproduction_level(params)),
                 tolerance = 1e-6)
})

test_that("getR0() does not need the species to be present", {
    params <- r0_params("centred", TRUE, TRUE)
    r0 <- getR0(params)
    initialN(params)[1, ] <- 0
    # The ratio of RDI to required recruitment is now 0 / 0
    expect_equal(getRDI(params), c("Target species" = 0))
    expect_equal(getR0(params), r0)

    params <- NS_params_small
    initialN(params)["Cod", ] <- 0
    r0 <- getR0(params, species = "Cod")
    expect_true(is.finite(r0) && r0 > 0)
})

test_that("getR0() uses the rates it is given", {
    params <- NS_params_small
    r0 <- getR0(params)
    expect_identical(getR0(params, e_growth = getEGrowth(params),
                           mort = getMort(params),
                           diffusion = getDiffusion(params),
                           e_repro = getERepro(params)),
                     r0)
    # R0 is proportional to the investment into reproduction ...
    expect_equal(getR0(params, e_repro = 2 * getERepro(params)), 2 * r0)
    # ... and falls as mortality rises
    expect_true(all(getR0(params, mort = 1.5 * getMort(params)) < r0))
})

test_that("getR0() selects species", {
    params <- NS_params_small
    r0 <- getR0(params)
    expect_named(r0, params@species_params$species)
    expect_identical(getR0(params, species = c("Cod", "Sprat")),
                     r0[c("Cod", "Sprat")])
    expect_identical(getR0(params, species = 2), r0[2])
})

test_that("getR0() and getNPerRecruit() check the shape of the rates", {
    params <- NS_params_small
    expect_error(getR0(params, mort = getMort(params)[, -1]), "`mort`")
    expect_error(getR0(params, e_repro = 1), "`e_repro`")
    expect_error(getNPerRecruit(params, e_growth = getEGrowth(params)[-1, ]),
                 "`e_growth`")
})

# getNPerRecruit ----
test_that("getNPerRecruit() is the upwind recursion in the default scheme", {
    params <- single_sp_params
    g <- getEGrowth(params)[1, ]
    mu <- getMort(params)[1, ]
    dw <- dw(params)
    i0 <- params@w_min_idx[[1]]
    # The upwind scheme carries abundance one size class past w_max
    top <- min(sum(w(params) <= species_params(params)$w_max) + 1,
               length(g))
    expected <- numeric(length(g))
    expected[i0] <- 1 / (g[i0] + mu[i0] * dw[i0])
    for (j in (i0 + 1):top) {
        expected[j] <- g[j - 1] * expected[j - 1] / (g[j] + mu[j] * dw[j])
    }
    expect_equal(getNPerRecruit(params)[1, ], expected, ignore_attr = TRUE)
})

test_that("getNPerRecruit() is the steady state per unit recruitment", {
    params <- r0_params("centred", TRUE, TRUE)
    n1 <- getNPerRecruit(params)
    expect_s3_class(n1, "ArraySpeciesBySize")
    expect_identical(dimnames(n1), dimnames(initialN(params)))
    expect_identical(attr(n1, "units"), "year/g")

    # The steady state is the abundance per recruit times the recruitment
    expect_equal(n1 * getRequiredRDD(params), initialN(params),
                 tolerance = 1e-10, ignore_attr = TRUE)
    # and R0 is the egg production of the abundance per recruit. This needs
    # the rates not to depend on the abundance, as here.
    expect_equal(getR0(params), getRDI(params, n = n1 * 1))
})
