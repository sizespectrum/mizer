test_that("We can set and get resource parameters", {
  params <- NS_params
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
  params <- setResource(NS_params, resource_rate = rr, 
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
  expect_error(setResource(NS_params, resource_dynamics = "fake"),
               'The resource dynamics function "fake" is not defined.')
    expect_error(setResource(NS_params, 
                             resource_capacity = 1, resource_level = 2),
                 "You should specify only either")
    expect_error(setResource(NS_params, 
                             resource_capacity = 1, resource_rate = 2),
                 "You should only provide either")
})



test_that("Can get and set resource_capacity slot", {
  params <- NS_params
  new <- 2 * resource_capacity(params)
  comment(new) <- "test"
  resource_capacity(params) <- new
  expect_identical(resource_capacity(params), new)
})
test_that("Can get and set resource_rate slot", {
  params <- NS_params
  new <- 2 * resource_rate(params)
  comment(new) <- "test"
  resource_rate(params) <- new
  expect_identical(resource_rate(params), new)
})
test_that("Can get and set resource_dynamics slot", {
    params <- NS_params
    new <- "resource_constant"
    comment(new) <- "test"
    resource_dynamics(params) <- new
    expect_identical(resource_dynamics(params), new)
})
