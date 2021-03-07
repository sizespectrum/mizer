test_that("compareParams", {
  params <- NS_params
  expect_equal(compareParams(params, params), character())
  # Change a species param
  params2 <- params
  params2@species_params$gamma[[1]] <- 1
  expect_equal(compareParams(params, params2)[[1]],
      'The following species parameters differ: Component "gamma": Mean absolute difference: 1')
  # Add a species param
  params@species_params$extra <- 1
  expect_true("params1 has the following additional species parameters: extra" %in%
                  compareParams(params, params2))
  params2@species_params$extra2 <- 2
  expect_true("params2 has the following additional species parameters: extra2" %in%
                  compareParams(params, params2))
  # Change a resource param
  params2@resource_params$lambda <- 1
  expect_true('The following resource parameters differ: Component "lambda": Mean absolute difference: 1.133333' %in%
                  compareParams(params, params2))
  # Change a slot
  params2@metab[] <- 1
  expect_true("The metab slots do not agree: Mean absolute difference: 896.7543" %in%
                  compareParams(params, params2))

  # Change size bins
  params2@w[[1]] <- 1
  expect_true("The community size bins differ." %in%
                  compareParams(params, params2))
  params2@w_full[[1]] <- 1
  expect_true("The resource size bins differ." %in%
                  compareParams(params, params2))
  # Add size bins
  params2@w <- c(1, params2@w)
  expect_true("The number of community size bins is different." %in%
                  compareParams(params, params2))
  params2@w_full <- c(1, params2@w_full)
  expect_true("The number of resource size bins is different." %in%
                  compareParams(params, params2))
})
