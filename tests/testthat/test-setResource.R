test_that("We can set and get resource parameters", {
  params <- NS_params
  rp <- resource_params(params)
  expect_identical(rp, params@resource_params)
  # Check that when called without parameters it leaves the params untouched
  # using expect_equal instead of expect_identical because the latter failed
  # on win_builder
  expect_equal(params, setResource(params))
  # Create example parameters
  resource_rate <- params@w_full
  comment(resource_rate) <- "resource_rate"
  names(resource_rate) <- names(params@rr_pp)
  resource_capacity <- 2 * params@w_full
  comment(resource_capacity) <- "resource_capacity"
  names(resource_capacity) <- names(params@cc_pp)
  r_pp <- 12
  comment(r_pp) <- "r_pp"
  kappa <- 13
  comment(kappa) <- "kappa"
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
                        r_pp = r_pp, kappa = kappa, lambda = lambda, n = n,
                        w_pp_cutoff = w_pp_cutoff, 
                        resource_dynamics = resource_dynamics)
  expect_identical(resource_params(params)$r_pp, r_pp)
  expect_identical(resource_params(params)$kappa, kappa)
  expect_identical(resource_params(params)$lambda, lambda)
  expect_identical(resource_params(params)$n, n)
  expect_identical(resource_params(params)$w_pp_cutoff, w_pp_cutoff)
  expect_identical(getResourceDynamics(params), resource_dynamics)
  expect_identical(getResourceRate(params), resource_rate)
  expect_identical(getResourceCapacity(params), resource_capacity)
  # Check that comments protect
  expect_message(params <- setResource(params, r_pp = r_pp, kappa = kappa,
                                       lambda = lambda, n = n,
                                       w_pp_cutoff = w_pp_cutoff),
                 "The resource carrying capacity has been commented")
  comment(params@cc_pp) <- NULL
  expect_message(params <- setResource(params, r_pp = r_pp, kappa = kappa,
                                       lambda = lambda, n = n,
                                       w_pp_cutoff = w_pp_cutoff),
                 "The resource intrinsic growth rate has been commented")
  # Check that setResource calculates rates correctly
  params <- setResource(NS_params, r_pp = r_pp, kappa = kappa,
                        lambda = lambda, n = n,
                        w_pp_cutoff = w_pp_cutoff)
  expect_identical(unname(getResourceRate(params)),
                   r_pp * params@w_full^(n - 1))
  resource_capacity <- getResourceCapacity(params)
  cc_pp <- kappa*params@w_full^(-lambda)
  expect_equivalent(cc_pp[params@w_full <= w_pp_cutoff],
                   resource_capacity[params@w_full <= w_pp_cutoff])
  # resource_params<- sets rates correctly
  resource_params(NS_params) <- resource_params(params)
  expect_identical(NS_params, params)
})
