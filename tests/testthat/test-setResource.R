test_that("We can set and get resource parameters", {
  params <- NS_params_small
  rp <- params@resource_params
  # Check that when called without parameters it leaves the params untouched
  # using expect_equal instead of expect_identical because the latter failed
  # on win_builder
  expect_unchanged(params, setResource(params))
  
  # Create example parameters
  resource_rate <- params@w_full
  comment(resource_rate) <- "resource_rate"
  names(resource_rate) <- names(params@rr_pp)
  resource_capacity <- 2 * params@w_full
  comment(resource_capacity) <- "resource_capacity"
  names(resource_capacity) <- names(params@cc_pp)
  rr <- 12
  comment(rr) <- "rr"
  cc <- 13
  comment(cc) <- "cc"
  lambda <- 1.79
  comment(lambda) <- "lambda"
  n <- 0.68
  comment(n) <- "n"
  w_pp_cutoff <- 0.5
  comment(w_pp_cutoff) <- "w_pp_cutoff"
  resource_dynamics <- "resource_constant"
  comment(resource_dynamics) <- "resource_dynamics"
  params <- setResource(params, resource_rate = resource_rate, 
                        resource_capacity = resource_capacity,
                        lambda = lambda, n = n,
                        w_pp_cutoff = w_pp_cutoff, 
                        resource_dynamics = resource_dynamics,
                        balance = FALSE)
  expect_identical(params@resource_params$lambda, lambda)
  expect_identical(params@resource_params$n, n)
  expect_identical(params@resource_params$w_pp_cutoff, w_pp_cutoff)
  expect_identical(params@resource_dynamics, resource_dynamics)
  expect_identical(params@rr_pp, resource_rate)
  expect_identical(params@cc_pp, resource_capacity)
  
  # Check that setResource calculates rates correctly
  params <- setResource(NS_params_small, resource_rate = rr, 
                        resource_capacity = cc,
                        lambda = lambda, n = n,
                        w_pp_cutoff = w_pp_cutoff,
                        balance = FALSE)
  expected <- rr * params@w_full^(n - 1)
  comment(expected) <- comment(rr)
  expect_identical(unname(params@rr_pp), expected)
  
  expected <- cc * params@w_full^(-lambda)
  comment(expected) <- comment(cc)
  expect_equal(expected[params@w_full <= w_pp_cutoff],
               unname(params@cc_pp[params@w_full <= w_pp_cutoff]))
})

test_that("setResource gives error", {
  expect_error(setResource(NS_params_small, resource_dynamics = "fake"),
               'The resource dynamics function "fake" is not defined.')
    expect_error(setResource(NS_params_small, 
                             resource_capacity = 1, resource_level = 2),
                 "You should specify only either")
    expect_error(setResource(NS_params_small, 
                             resource_capacity = 1, resource_rate = 2),
                 "You should only provide either")
})



test_that("Can get and set resource_capacity slot", {
  params <- NS_params_small
  new <- 2 * resource_capacity(params)
  comment(new) <- "test"
  resource_capacity(params) <- new
  expect_equal(resource_capacity(params), new, ignore_attr = TRUE)
  expect_identical(comment(resource_capacity(params)), "test")
})
test_that("Can get and set resource_rate slot", {
  params <- NS_params_small
  new <- 2 * resource_rate(params)
  comment(new) <- "test"
  resource_rate(params) <- new
  expect_equal(resource_rate(params), new, ignore_attr = TRUE)
  expect_identical(comment(resource_rate(params)), "test")
})
test_that("resource accessors return documented values", {
    params <- NS_params_small
    expect_equal(resource_rate(params), params@rr_pp, ignore_attr = TRUE)
    expect_equal(resource_capacity(params), params@cc_pp, ignore_attr = TRUE)
    expect_equal(resource_level(params), params@initial_n_pp / params@cc_pp,
                 ignore_attr = TRUE)
})
test_that("Can get and set resource_level slot", {
    params <- NS_params_small
    new <- rep(0.5, length(resource_level(params)))
    resource_level(params) <- new
    expected_level <- new
    expected_level[initialNResource(params) == 0] <- NaN
    expected_capacity <- initialNResource(params) / expected_level
    expected_capacity[is.nan(expected_capacity)] <- 0
    expect_equal(resource_level(params), expected_level, ignore_attr = TRUE)
    expect_equal(resource_capacity(params), expected_capacity,
                 ignore_attr = TRUE)
})

test_that("resource_level can be set to 1", {
    params <- NS_params_small
    expect_warning(
        params_new <- setResource(params, resource_level = 1),
        "division by zero"
    )

    expect_true(all(is.finite(resource_rate(params_new))))
    expect_true(all(resource_capacity(params_new) >= initialNResource(params)))
})

test_that("Can get and set resource_dynamics slot", {
    params <- NS_params_small
    new <- "resource_constant"
    comment(new) <- "test"
    resource_dynamics(params) <- new
    expect_identical(resource_dynamics(params), new)
})

test_that("w_pp_cutoff can be decreased without providing carrying_capacity", {
    params <- NS_params_small
    old_cutoff <- resource_params(params)$w_pp_cutoff
    # Choose a new cutoff that is smaller
    new_cutoff <- old_cutoff / 2
    
    # This should work and cut off both carrying capacity and initial abundance
    params_new <- setResource(params, w_pp_cutoff = new_cutoff)
    
    # Check that w_pp_cutoff was updated
    expect_equal(resource_params(params_new)$w_pp_cutoff, new_cutoff)
    
    # Check that carrying capacity is zero above the new cutoff
    w_full <- params_new@w_full
    expect_true(all(params_new@cc_pp[w_full >= new_cutoff] == 0))
    
    # Check that initial resource abundance is zero above the new cutoff
    expect_true(all(params_new@initial_n_pp[w_full >= new_cutoff] == 0))
    
    # Check that carrying capacity is non-zero below the new cutoff
    # (at least somewhere)
    expect_true(any(params_new@cc_pp[w_full < new_cutoff] > 0))
})

test_that("w_pp_cutoff cannot be increased without providing carrying_capacity", {
    params <- NS_params_small
    old_cutoff <- resource_params(params)$w_pp_cutoff
    # Try to increase the cutoff
    new_cutoff <- old_cutoff * 2
    
    # This should give an error
    expect_error(
        setResource(params, w_pp_cutoff = new_cutoff),
        "You cannot increase w_pp_cutoff without also providing the resource_capacity"
    )
})

test_that("w_pp_cutoff can be changed when providing carrying_capacity", {
    params <- NS_params_small
    old_cutoff <- resource_params(params)$w_pp_cutoff
    new_cutoff <- old_cutoff * 2
    
    # This should work when we also provide carrying capacity
    params_new <- setResource(params, w_pp_cutoff = new_cutoff, 
                             resource_capacity = 10, balance = FALSE)
    
    # Check that w_pp_cutoff was updated
    expect_equal(resource_params(params_new)$w_pp_cutoff, new_cutoff)
    
    # Check that carrying capacity follows the new cutoff
    w_full <- params_new@w_full
    expect_true(all(params_new@cc_pp[w_full >= new_cutoff] == 0))
})

test_that("second_order_w bin-averages the resource rate and capacity", {
    rr <- 12
    cc <- 13
    lambda <- 2.05
    n <- 0.68
    w_pp_cutoff <- 0.5

    base <- setResource(NS_params_small, resource_rate = rr,
                        resource_capacity = cc, lambda = lambda, n = n,
                        w_pp_cutoff = w_pp_cutoff, balance = FALSE)
    p2 <- base
    second_order_w(p2) <- c(bin_average = TRUE)
    p2 <- setResource(p2, resource_rate = rr, resource_capacity = cc,
                      lambda = lambda, n = n, w_pp_cutoff = w_pp_cutoff,
                      balance = FALSE)

    w <- p2@w_full
    dw <- p2@dw_full

    # Rate: exact bin average of rr * w^(n-1).
    rate_expected <- rr * power_law_bin_average(w, dw, n - 1)
    expect_equal(unname(p2@rr_pp), rate_expected)
    # Capacity: exact bin average of cc * w^(-lambda), cut at w_pp_cutoff.
    cap_expected <- cc * power_law_bin_average(w, dw, -lambda,
                                               w_max = w_pp_cutoff)
    expect_equal(unname(p2@cc_pp), cap_expected)

    # The bin averages differ from the left-edge point values (the power law
    # varies across the bin).
    expect_false(isTRUE(all.equal(unname(p2@rr_pp), unname(base@rr_pp))))
    expect_false(isTRUE(all.equal(unname(p2@cc_pp), unname(base@cc_pp))))
})

test_that("second_order_w capacity straddling-bin gets the partial average", {
    cc <- 13
    lambda <- 2.05
    w_pp_cutoff <- 0.5
    p2 <- NS_params_small
    second_order_w(p2) <- c(bin_average = TRUE)
    p2 <- setResource(p2, resource_capacity = cc, lambda = lambda,
                      w_pp_cutoff = w_pp_cutoff, balance = FALSE)
    w <- p2@w_full
    dw <- p2@dw_full

    # Bins entirely above the cutoff are zero.
    expect_true(all(p2@cc_pp[w >= w_pp_cutoff] == 0))
    # The straddling bin (cutoff strictly inside) gets the partial power-law
    # average over the part below the cutoff, not the full-bin average.
    straddle <- which(w < w_pp_cutoff & (w + dw) > w_pp_cutoff)
    if (length(straddle) == 1) {
        partial <- cc * (w_pp_cutoff^(1 - lambda) - w[straddle]^(1 - lambda)) /
            ((1 - lambda) * dw[straddle])
        expect_equal(unname(p2@cc_pp[straddle]), partial)
        full <- cc * power_law_bin_average(w[straddle], dw[straddle], -lambda)
        expect_lt(unname(p2@cc_pp[straddle]), full)
    }
})

test_that("default resource construction is byte-identical (first order)", {
    # With second_order_w off (default) the rate and capacity are the plain
    # left-edge power laws, unchanged from previous mizer.
    rr <- 12
    cc <- 13
    lambda <- 2.05
    n <- 0.68
    w_pp_cutoff <- 0.5
    p <- setResource(NS_params_small, resource_rate = rr,
                     resource_capacity = cc, lambda = lambda, n = n,
                     w_pp_cutoff = w_pp_cutoff, balance = FALSE)
    w <- p@w_full
    expect_equal(unname(p@rr_pp), rr * w^(n - 1))
    cap <- cc * w^(-lambda)
    cap[w >= w_pp_cutoff] <- 0
    expect_equal(unname(p@cc_pp), cap)
})

test_that("resource_params<- rebuilds capacity and rate arrays (issue #439)", {
    p <- NS_params_small
    cc0 <- resource_capacity(p)
    rr0 <- resource_rate(p)

    # 1. resource_params(p)$kappa <- x scales cc_pp
    p1 <- p
    resource_params(p1)$kappa <- 10 * resource_params(p1)$kappa
    expect_equal(as.numeric(resource_capacity(p1)), 10 * as.numeric(cc0))
    # resource_params<- does not balance, so rr_pp is left unchanged
    expect_equal(as.numeric(resource_rate(p1)), as.numeric(rr0))

    # 2. resource_params(p)$lambda <- x changes the slope of cc_pp
    p2 <- p
    rp2 <- resource_params(p2)
    rp2$kappa <- 1e14
    rp2$lambda <- 2.2
    resource_params(p2) <- rp2
    expect_equal(resource_params(p2)$lambda, 2.2)
    expect_false(isTRUE(all.equal(as.numeric(resource_capacity(p2)), as.numeric(cc0))))

    # 3. Freeze mechanism: manually set resource_capacity is preserved when scalar changes
    p3 <- p
    custom_cc <- 2 * resource_capacity(p3)
    resource_capacity(p3) <- custom_cc
    expect_identical(comment(p3@cc_pp), "set manually")
    
    # Modifying kappa now does not overwrite custom cc_pp array and, because
    # resource_params<- no longer balances, the freeze comment survives too.
    resource_params(p3)$kappa <- 10 * resource_params(p3)$kappa
    expect_equal(resource_capacity(p3), custom_cc, ignore_attr = TRUE)
    expect_identical(comment(p3@cc_pp), "set manually")
    
    # setResource(..., reset = TRUE) resets to scalar-driven capacity
    p4 <- setResource(p3, reset = TRUE)
    expect_null(comment(p4@cc_pp))
    expect_equal(as.numeric(resource_capacity(p4)), 10 * as.numeric(cc0))

    # 4. w_pp_cutoff change via resource_params<- truncates capacity
    p5 <- p
    old_cutoff <- resource_params(p5)$w_pp_cutoff
    new_cutoff <- old_cutoff / 2
    resource_params(p5)$w_pp_cutoff <- new_cutoff
    expect_equal(resource_params(p5)$w_pp_cutoff, new_cutoff)
    w_full <- p5@w_full
    expect_true(all(p5@cc_pp[w_full >= new_cutoff] == 0))
    expect_true(all(p5@initial_n_pp[w_full >= new_cutoff] == 0))
})

test_that("resource_params<- does not balance and rate-side scalars take effect", {
    p <- NS_params_small
    # Give the model a scalar r_pp so the rate is driven by resource_params
    p <- setResource(p, resource_rate = 1, balance = FALSE)
    rr1 <- as.numeric(resource_rate(p))

    # Changing the rate-side scalar r_pp now updates rr_pp (was inert before
    # because balancing re-derived the rate from the capacity)
    resource_params(p)$r_pp <- 2
    expect_equal(as.numeric(resource_rate(p)), 2 * rr1)

    # Sequential scalar edits accumulate: a later capacity-side change does not
    # discard the earlier rate-side change (mirrors species_params<-)
    cc_before <- as.numeric(resource_capacity(p))
    rr_before <- as.numeric(resource_rate(p))
    resource_params(p)$kappa <- 5 * resource_params(p)$kappa
    expect_equal(as.numeric(resource_rate(p)), rr_before)
    expect_equal(as.numeric(resource_capacity(p)), 5 * cc_before)
})

test_that("resource setters take a balance argument", {
    p <- NS_params_small
    rr0 <- as.numeric(resource_rate(p))

    # balance = FALSE leaves the opposite side untouched
    p_f <- p
    resource_capacity(p_f, balance = FALSE) <- 2 * resource_capacity(p_f)
    expect_equal(as.numeric(resource_rate(p_f)), rr0)

    # the default (balance = TRUE) rebalances the rate to the new capacity
    p_t <- p
    resource_capacity(p_t) <- 2 * resource_capacity(p_t)
    expect_false(isTRUE(all.equal(as.numeric(resource_rate(p_t)), rr0)))
})

test_that("incidental balancing does not overwrite a frozen array (freeze wins)", {
    p <- NS_params_small
    # Freeze the capacity manually
    resource_capacity(p) <- 2 * resource_capacity(p)
    expect_false(is.null(comment(p@cc_pp)))
    cc_frozen <- as.numeric(resource_capacity(p))

    # An incidental balance (no rate/capacity supplied) would re-derive the
    # capacity, but the frozen capacity must be preserved and a warning issued.
    expect_warning(
        p2 <- setResource(p, balance = TRUE),
        "set manually"
    )
    expect_equal(as.numeric(resource_capacity(p2)), cc_frozen)
    expect_false(is.null(comment(p2@cc_pp)))
})

test_that("explicit balancing overrides a frozen complementary array", {
    p <- NS_params_small
    # Freeze the rate manually
    resource_rate(p) <- resource_rate(p)
    expect_false(is.null(comment(p@rr_pp)))

    # Explicitly supplying a new capacity with balance = TRUE (the default)
    # re-derives the rate, overriding the freeze, so the resource stays at
    # steady state.
    p2 <- setResource(p, resource_capacity = 2 * resource_capacity(p))
    N <- initialNResource(p2)
    mu <- getResourceMort(p2)
    r <- p2@rr_pp
    c <- p2@cc_pp
    sel <- (mu + r) > 0
    expect_equal((r * c / (mu + r))[sel], N[sel], ignore_attr = TRUE)

    # ... but balance = FALSE keeps the frozen rate untouched.
    p3 <- setResource(p, resource_capacity = 2 * resource_capacity(p),
                      balance = FALSE)
    expect_equal(resource_rate(p3), resource_rate(p), ignore_attr = TRUE)
})
