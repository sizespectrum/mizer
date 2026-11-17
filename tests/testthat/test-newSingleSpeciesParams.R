# newSingleSpeciesParams ----
test_that("newSingleSpeciesParams works", {
    params <- newSingleSpeciesParams()
    no_w <- length(params@w)
    expect_equal(dim(params@initial_n), c(1, no_w))
    expect_equal(params@w[1], params@species_params$w_min, ignore_attr = TRUE)
    expect_equal(params@w[no_w], params@species_params$w_max, ignore_attr = TRUE)
    sim <- project(params, t_max = 1)
    expect_equal(sim@n[1, 1, ], sim@n[2, 1, ])
    expect_error(newSingleSpeciesParams(f0 = 1, gamma = 1, info_level = 0),
                 "`f0` must be finite and in the interval")
})

test_that("newSingleSpeciesParams uses the scheme of the defaults edition", {
    withr::local_options(mizer_defaults_edition = 2)
    expect_identical(second_order_w(newSingleSpeciesParams(info_level = 0)),
                     list(flux = "upwind", bin_average = FALSE))
    withr::local_options(mizer_defaults_edition = 3)
    params <- newSingleSpeciesParams(info_level = 0)
    expect_identical(second_order_w(params),
                     list(flux = "van_leer", bin_average = TRUE))
    p_true <- newSingleSpeciesParams(second_order_w = TRUE, info_level = 0)
    p_true@time_created <- params@time_created
    p_true@time_modified <- params@time_modified
    expect_equal(params, p_true)
    expect_identical(
        second_order_w(newSingleSpeciesParams(second_order_w = FALSE,
                                              info_level = 0)),
        list(flux = "upwind", bin_average = FALSE))
})

test_that("newSingleSpeciesParams documents and applies grid and deprecation behaviour", {
    expect_message(
        params <- newSingleSpeciesParams(w_max = 100, w_min = 0.001, no_w = 2),
        "Increased no_w to 26"
    )
    expect_equal(length(w(params)), 26)

    expect_snapshot_warning(params2 <- newSingleSpeciesParams(R_factor = Inf))
    expect_equal(unname(reproduction_level(params2)), rep(0, nrow(species_params(params2))))
})

test_that("Sets given_species_params", {
    # Calling `given_species_params<-()` should not make a change
    params <- newSingleSpeciesParams()
    p2 <- params
    given_species_params(p2) <- given_species_params(p2)
    expect_unchanged(p2, params)
})

test_that("newSingleSpeciesParams includes external diffusion", {
    for (sow in c(FALSE, TRUE)) {
        params <- newSingleSpeciesParams(D_ext = 0.1, second_order_w = sow,
                                         info_level = 0)
        expect_equal(unname(species_params(params)$D_ext), 0.1)
        expect_equal(unname(given_species_params(params)$D_ext), 0.1)
        exponent <- species_params(params)$n[[1]] + 1
        expected <- if (sow) {
            ((params@w + params@dw)^(exponent + 1) -
                 params@w^(exponent + 1)) / ((exponent + 1) * params@dw)
        } else {
            params@w^exponent
        }
        expect_equal(as.numeric(params@ext_diffusion), 0.1 * expected)
        expect_null(comment(params@ext_diffusion))
        sim <- project(params, t_max = 0.01, dt = 0.001, t_save = 0.01,
                       effort = 0, progress_bar = FALSE)
        expect_true(all(is.finite(sim@n)))
        expect_true(all(sim@n >= 0))
        if (!sow) {
            # The first-order initial spectrum must balance diffusion as well
            # as growth and mortality, including at the egg boundary.
            expect_equal(sim@n[2, 1, ] / sim@n[1, 1, ],
                         rep(1, length(params@w)), ignore_attr = TRUE,
                         tolerance = 1e-10)
        }
    }
})

test_that("newSingleSpeciesParams defaults to zero external diffusion", {
    # single_sp_params is built without a D_ext argument.
    expect_true(all(single_sp_params@ext_diffusion == 0))
    expect_equal(unname(given_species_params(single_sp_params)$D_ext), 0)
    params <- newSingleSpeciesParams(D_ext = 0, info_level = 0)
    expect_equal(params@initial_n, single_sp_params@initial_n)
    expect_equal(params@species_params, single_sp_params@species_params)
})

test_that("newSingleSpeciesParams keeps initial_n inside the projection support", {
    # With w_max = 77.7 the top grid point rounds to just above w_max, so the
    # second-order support stops one bin below the upwind support used for
    # construction.
    for (D_ext in c(0, 0.1)) {
        params <- newSingleSpeciesParams(w_max = 77.7, D_ext = D_ext,
                                         second_order_w = TRUE,
                                         info_level = 0)
        top <- support_top_idx(params)
        expect_lt(top, length(params@w))
        expect_true(params@initial_n[1, top] > 0)
        expect_true(all(params@initial_n[1, -seq_len(top)] == 0))
    }
})

test_that("single-species diffusive equilibrium also uses bin-averaged rates", {
    params <- newSingleSpeciesParams(D_ext = 0.1, info_level = 0,
                                     second_order_w = list(flux = "upwind",
                                                           bin_average = TRUE))
    sim <- project(params, t_max = 0.01, dt = 0.001, t_save = 0.01,
                   effort = 0, progress_bar = FALSE)
    expect_equal(sim@n[2, 1, ] / sim@n[1, 1, ], rep(1, length(params@w)),
                 ignore_attr = TRUE, tolerance = 1e-10)
})

test_that("newSingleSpeciesParams validates the diffusion coefficient", {
    for (D_ext in list(-1, NA_real_, Inf, NaN, c(0, 0.1), "0.1", NULL)) {
        expect_error(newSingleSpeciesParams(D_ext = D_ext, info_level = 0),
                     "D_ext")
    }
})
