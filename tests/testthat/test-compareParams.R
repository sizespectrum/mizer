test_that("compareParams", {
  local_reproducible_output()
  params <- NS_params
  expect_equal(compareParams(params, params), "No differences")
  # Change a species param
  params2 <- params
  params2@species_params$gamma[[1]] <- 1
  expect_equal(compareParams(params, params2)[[1]],
      'The following species parameters differ: Component "gamma": Mean absolute difference: 1')
  # Add a species param
  params1 <- params
  params1@species_params$extra <- 1
  expect_true("params1 has the following additional species parameters: extra" %in%
                  compareParams(params1, params))
  params2@species_params$extra2 <- 2
  expect_true("params2 has the following additional species parameters: extra2" %in%
                  compareParams(params, params2))
  # Change a given species param
  params2 <- params
  params2@given_species_params$gamma[[1]] <- 1
  expect_equal(compareParams(params, params2)[[1]],
               'The following given species parameters differ: Component "gamma": Mean absolute difference: 1')
  # Add a given species param
  params1 <- params
  params1@given_species_params$extra <- 1
  expect_true("params1 has the following additional given species parameters: extra" %in%
                  compareParams(params1, params))
  params2@given_species_params$extra2 <- 2
  expect_true("params2 has the following additional given species parameters: extra2" %in%
                  compareParams(params, params2))
  # Change a resource param
  params2@resource_params$lambda <- 1
  expect_true('The following resource parameters differ: Component "lambda": Mean absolute difference: 1.133333' %in%
                  compareParams(params, params2))
  # Change a slot
  params2@metab[] <- 1
  expect_true("The metab slots do not agree: Mean absolute difference: 896.7543" %in%
                  compareParams(params, params2))

  # Change size bin at small sizes
  params2@w_full[[1]] <- 1
  expect_true("The resource size bins differ." %in%
                  compareParams(params, params2))
  # Change a size bin in the fish range
  fish_offset <- length(params2@w_full) - length(params2@w)
  params2@w[[10]] <- 1
  params2@w_full[[10 + fish_offset]] <- 1
  expect_true("The community size bins differ." %in%
                  compareParams(params, params2))
})
