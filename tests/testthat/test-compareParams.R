test_that("compareParams", {
  local_reproducible_output()
  sink(nullfile())
  on.exit(sink(), add = TRUE, after = FALSE)
  params <- NS_params_small
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
  expect_true('The following resource parameters differ: Component "lambda": Mean absolute difference: 1.05' %in%
                  compareParams(params, params2))
  # Change a slot
  params2@metab[] <- 1
  expect_true(any(startsWith(compareParams(params, params2),
                             "The metab slots do not agree: Mean absolute difference: 687.7477")))

  # Change only the comment attribute of a slot
  params2 <- params
  comment(params2@metab) <- "a comment"
  expect_true(grepl("comment", compareParams(params, params2)[[1]],
                    ignore.case = TRUE))

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

test_that("compareParams reports mismatched species and gear counts", {
  local_reproducible_output()
  sink(nullfile())
  on.exit(sink(), add = TRUE, after = FALSE)
  params <- NS_params_small

  # Different number of species: report the count difference and the extra
  # species, but do not try to compare the incompatible arrays.
  params_fewer <- removeSpecies(params, species_params(params)$species[[1]])
  result <- compareParams(params, params_fewer)
  expect_true("The number of species is different." %in% result)
  expect_true(any(startsWith(result, "params1 has the following additional species:")))
  expect_false(any(startsWith(result, "The metab slots do not agree")))
  # Identical messages are not repeated
  expect_equal(anyDuplicated(result), 0L)

  # Different number of gears
  params_gears <- params
  gear_params(params_gears)$gear <- "single"
  expect_true("The number of gears is different." %in%
                  compareParams(params, params_gears))
})
