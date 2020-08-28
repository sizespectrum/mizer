params <- NS_params

# validGearParams ----
test_that("validGearParams works", {
    sp <- validSpeciesParams(
        data.frame(species = c("species1", "species2"),
                   w_inf = c(100, 1000),
                   stringsAsFactors = FALSE))
    # gear_params is allowed to have zero rows
    gp <- validGearParams(data.frame(), sp)
    expect_identical(gp, data.frame(species = list(), gear = list()))
    # There must be columns `species` and `gear`
    gp <- data.frame(species = 1)
    expect_error(validGearParams(gp, sp), 
                 "`gear_params` must have columns 'species' and 'gear'.")
    # Any species-gear pair is allowed to appear at most once
    gp <- data.frame(species = c("species1", "species1"), gear = c("g", "g"),
                     stringsAsFactors = FALSE)
    expect_error(validGearParams(gp, sp), 
                 "Some species - gear pairs appear more than once.")
    # Any species that appears must also appear in the `species_params` data frame.
    gp <- data.frame(species = c("species1", "species3"), gear = c("g", "g"),
                     stringsAsFactors = FALSE)
    expect_error(validGearParams(gp, sp), 
                 "The gear_params dataframe contains species that do not exist in the model.")
    # There must be a `sel_fun` column
    gp <- validGearParams(
        data.frame(species = c("species1", "species2"), gear = c("g", "g"),
                   stringsAsFactors = FALSE),
        sp)
    expect_identical(gp$sel_func, c("knife_edge", "knife_edge"))
    expect_identical(gp$knife_edge_size, c(25, 250))
    # There must be a catchability column
    expect_identical(gp$catchability, c(1, 1))
    # Defaults for NAs
    gp$gear[[1]] <- NA
    expect_identical(validGearParams(gp, sp)$gear[[1]], "species1")
    gp$sel_func[[1]] <- NA
    expect_identical(validGearParams(gp, sp)$sel_func[[1]], "knife_edge")
    gp$catchability[[2]] <- NA
    expect_identical(validGearParams(gp, sp)$catchability[[2]], 1)
    gp$knife_edge_size[[2]] <- NA
    expect_identical(validGearParams(gp, sp)$knife_edge_size[[2]], 250)
})

# validEffortVector ----
test_that("validEffort works", {
    params <- NS_params
    ie <- params@initial_effort
    # an already valid vector is not changed
    expect_identical(validEffortVector(ie, params), ie)
    # A scrambled vector is put in the right order
    ies <- ie[c(2,3,1,4)]
    expect_identical(validEffortVector(ies, params), ie)
    # A single number is converted into a constant vector
    ie[] <- 2
    expect_identical(validEffortVector(2, params), ie)
    # The length of the vector is checked
    expect_error(validEffortVector(ie[1:3], params),
                 "Effort vector must be the same length as the number of fishing gears.")
    # The names are checked
    names(ie)[[1]] <- "test"
    expect_error(validEffortVector(ie, params), 
                 "Gear names in the MizerParams object")
})
test_that("validEffortParams works when no gears are set up", {
    params <- newMultispeciesParams(NS_species_params,
                                    gear_params = data.frame())
    expect_length(validEffortVector(1, params), 0)
    expect_length(validEffortVector(NULL, params), 0)
})

# setFishing and gear_params ----
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

test_that("Duplicate gear-species pairs give error", {
    gp <- rbind(NS_params@gear_params, NS_params@gear_params[1, ])
    expect_error(gear_params(params) <- gp,
                 "Some species - gear pairs appear more than once")
})

test_that("Non-existing species give error", {
    gp <- NS_params@gear_params
    gp$species <- as.character(gp$species)
    gp$species[[1]] <- "test"
    expect_error(gear_params(params) <- gp,
                 "The gear_params dataframe contains species that do not exist in the model.")
})
