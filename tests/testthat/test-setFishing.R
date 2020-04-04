params <- NS_params

test_that("Set Fishing works", {
    expect_identical(gear_params(params), params@gear_params)
    expect_identical(params, setFishing(params))
    gear_params(params) <- params@gear_params
    expect_identical(params, NS_params)
})

test_that("Setting selectivity works", {
    selectivity <- getSelectivity(params)
    expect_identical(selectivity, params@selectivity)
    selectivity[1, 1, 1] <- 111
    comment(selectivity) <- "selectivity"
    params <- setFishing(params, selectivity = selectivity)
    expect_identical(params@selectivity[1, 1, 1], 111)
    expect_identical(comment(getSelectivity(params)), "selectivity")
})

test_that("Setting catchability works", {
    catchability <- getCatchability(params)
    expect_identical(catchability, params@catchability)
    catchability[2, 2] <- 22
    comment(catchability) <- "catchability"
    params <- setFishing(params, catchability = catchability)
    expect_identical(params@catchability[2, 2], 22)
    expect_identical(comment(getCatchability(params)), "catchability")
})

test_that("Comments protect slots", {
    catchability <- getCatchability(params)
    comment(catchability) <- "catchability"
    params <- setFishing(params, catchability = catchability)
    expect_message(gear_params(params) <- params@gear_params,
                   "The catchability has been commented")
    
    selectivity <- getSelectivity(params)
    comment(selectivity) <- "selectivity"
    expect_message(params <- setFishing(params, selectivity = selectivity),
                   "The catchability has been commented")
    expect_message(gear_params(params) <- params@gear_params,
                   "The selectivity has been commented")
})

test_that("Arguments of wrong dimension throw errors", {
    catchability <- getCatchability(params)
    expect_error(setFishing(params, catchability = catchability[1:3, ]))
                # "dim(catchability)[[1]] not equal to dim(selectivity)[[1]]"
    expect_error(setFishing(params, catchability = catchability[, 1:3]))
                 #"dim(catchability)[[2]] not equal to no_sp")
    
    selectivity <- getSelectivity(params)
    expect_error(setFishing(params, selectivity = selectivity[1:3, ]),
                 "incorrect number of dimensions")
    expect_error(setFishing(params, selectivity = selectivity[1:3, , ]))
                 #"dim(selectivity)[[1]] not equal to no_gears")
})

test_that("We get errors when gear parameters are missing", {
    gear_params <- params@gear_params
    expect_error(gear_params(params)$knife_edge_size[[2]] <- NA,
                 "Some selectivity parameters are NA")
    expect_error(gear_params(params)$knife_edge_size <- NULL,
                 "Some arguments needed for the selectivity function are missing")
})

test_that("Dimensions after number of gears has increased", {
    params <- NS_params
    gear_params(params)$gear <- params@species_params$species
    no_gears <- nrow(params@species_params)
    expect_identical(dim(params@catchability)[[1]], no_gears)
    expect_identical(dim(params@selectivity)[[1]], no_gears)
    expect_identical(dimnames(params@catchability)[[1]],
                     params@species_params$species)
    expect_identical(dimnames(params@selectivity)[[1]],
                     params@species_params$species)
})
