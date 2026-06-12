params <- NS_params_small

# validGearParams ----
test_that("validGearParams works", {
    sp <- validSpeciesParams(
        data.frame(species = c("species1", "species2"),
                   w_max = c(100, 1000),
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
    expect_identical(gp$catchability,
                     rep(ifelse(defaults_edition() < 2, 1, 0.3), 2))
    # Defaults for NAs
    gp$gear[[1]] <- NA
    expect_identical(validGearParams(gp, sp)$gear[[1]], "species1")
    gp$sel_func[[1]] <- NA
    expect_identical(validGearParams(gp, sp)$sel_func[[1]], "knife_edge")
    gp$catchability[[2]] <- NA
    expect_identical(validGearParams(gp, sp)$catchability[[2]],
                     ifelse(defaults_edition() < 2, 1, 0.3))
    gp$knife_edge_size[[2]] <- NA
    expect_identical(validGearParams(gp, sp)$knife_edge_size[[2]], 250)
    
    # The rownames must be of the form "species, gear"
    gp$species <- c("species1", "species1")
    gp$gear <- c("g1", "g2")
    expect_identical(rownames(validGearParams(gp, sp)), 
                     c("species1, g1", "species1, g2"))
    
    # Factors are converted to strings
    gp <- data.frame(species = factor("species1"), gear = factor("g"),
                     stringsAsFactors = TRUE)
    sp <- data.frame(species = "species1", w_max = 100)
    gp <- validGearParams(gp, sp)
    expect_identical(gp$species, "species1")
    expect_identical(gp$gear, "g")
})

# validEffortVector ----
test_that("validEffort works", {
    params <- NS_params_small
    ie <- params@initial_effort
    # an already valid vector is not changed
    expect_identical(validEffortVector(ie, params), ie)
    # A scrambled vector is put in the right order
    ies <- ie[c(2,3,1)]
    expect_identical(validEffortVector(ies, params), ie)
    # NA's are replaced by the edition-specific default
    ie[2] <- ifelse(defaults_edition() < 2, 0, 1)
    iesn <- ie
    iesn[2] <- NA
    expect_identical(validEffortVector(iesn, params), ie)
    # A single number is converted into a constant vector
    ie[] <- 2
    expect_identical(validEffortVector(2, params), ie)
    # A shortened vector is expanded with the edition-specific default
    expect_identical(validEffortVector(ie[c(1,2)], params)[[3]],
                     ifelse(defaults_edition() < 2, 0, 1))
                
    # The names are checked
    names(ie)[[1]] <- "test"
    expect_error(validEffortVector(ie, params), 
                 "it has names that are not among the gear names")
})
test_that("validEffortParams works when no gears are set up", {
    params <- newMultispeciesParams(NS_species_params_small,
                                    gear_params = data.frame(), info_level = 0)
    expect_length(validEffortVector(1, params), 0)
    expect_length(validEffortVector(NULL, params), 0)
})

# setFishing and gear_params ----
test_that("Set Fishing works", {
    params1 <- params
    expect_identical(gear_params(params), params@gear_params)
    expect_unchanged(params, setFishing(params))
    gear_params(params) <- params@gear_params
    expect_unchanged(params, params1)
    # has updated time_modified
    expect_false(identical(params@time_modified, params1@time_modified))
})

test_that("Setting selectivity works", {
    selectivity <- getSelectivity(params)
    expect_identical(selectivity, params@selectivity)
    expect_identical(selectivity(params), selectivity)
    selectivity[1, 1, 1] <- 111
    comment(selectivity) <- "selectivity"
    params <- setFishing(params, selectivity = selectivity)
    expect_identical(params@selectivity[1, 1, 1], 111)
    expect_identical(comment(getSelectivity(params)), "selectivity")
})

test_that("Setting catchability works", {
    catchability <- getCatchability(params)
    expect_identical(catchability, params@catchability)
    expect_identical(catchability(params), catchability)
    catchability[2, 2] <- 22
    comment(catchability) <- "catchability"
    params <- setFishing(params, catchability = catchability)
    expect_identical(params@catchability[2, 2], 22)
    expect_identical(comment(getCatchability(params)), "catchability")
})

test_that("Comment works on selectivity", {
    params <- NS_params_small
    # if no comment, it is set automatically
    selectivity <- params@selectivity
    params <- setFishing(params, selectivity = selectivity)
    expect_identical(comment(params@selectivity), "set manually")
    
    # comment is stored
    comment(selectivity) <- "test"
    params <- setFishing(params, selectivity = selectivity)
    expect_identical(comment(params@selectivity), "test")
    
    # if no comment, previous comment is kept
    comment(selectivity) <- NULL
    params <- setFishing(params, selectivity = selectivity)
    expect_identical(comment(params@selectivity), "test")
    
    # no message when nothing changes
    expect_message(setFishing(params), NA)
    # but message when a change is not stored due to comment
    params@gear_params$knife_edge_size <- 0
    expect_message(setFishing(params),  "has been commented")
    # Can reset
    p <- setFishing(params, reset = TRUE)
    expect_equal(p@selectivity[1, 1, 1], 1, ignore_attr = TRUE)
    expect_warning(setFishing(params, selectivity = selectivity,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

test_that("Comment works on catchability", {
    params <- NS_params_small
    # if no comment, it is set automatically
    catchability <- params@catchability
    params <- setFishing(params, catchability = catchability)
    expect_identical(comment(params@catchability), "set manually")
    
    # comment is stored
    comment(catchability) <- "test"
    params <- setFishing(params, catchability = catchability)
    expect_identical(comment(params@catchability), "test")
    
    # if no comment, previous comment is kept
    comment(catchability) <- NULL
    params <- setFishing(params, catchability = catchability)
    expect_identical(comment(params@catchability), "test")
    
    # no message when nothing changes
    expect_message(setFishing(params), NA)
    # but message when a change is not stored due to comment
    params@gear_params$catchability <- 2
    expect_message(setFishing(params),  "has been commented")
    # Can reset
    p <- setFishing(params, reset = TRUE)
    expect_equal(p@catchability[1, 1], 2, ignore_attr = TRUE)
    expect_warning(setFishing(params, catchability = catchability,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

test_that("We can change gears via catchability and selectivity arrays", {
    catchability <- getCatchability(params)
    sel <- c(1, 2)
    expect_error(setFishing(params, catchability = catchability[sel, ]),
                 "you also need to supply a selectivity array")
    selectivity <- getSelectivity(params)
    p2 <- setFishing(params, catchability = catchability[sel, ],
                     selectivity = selectivity[sel, , ])
    comment(p2@selectivity) <- NULL
    expect_identical(p2@selectivity, selectivity[sel, , ])
    comment(p2@catchability) <- NULL
    expect_identical(p2@catchability, catchability[sel, ])
    expect_identical(p2@initial_effort, params@initial_effort[sel])
    expect_error(setFishing(params, catchability = catchability[sel, ],
                            selectivity = selectivity),
                 "not equal to no_gears")
})

test_that("Arguments of wrong dimension throw errors", {
    catchability <- getCatchability(params)
    expect_error(setFishing(params, catchability = catchability[, 1:2]),
                 "not equal to no_sp")

    selectivity <- getSelectivity(params)
    expect_error(setFishing(params, selectivity = selectivity[1:2, ]),
                 "incorrect number of dimensions")
    expect_error(setFishing(params, selectivity = selectivity[1:2, , ]),
                 "not equal to no_gears")
})

test_that("Wrong dimnames throw errors or get fixed", {
    catchability <- getCatchability(params)
    cw <- catchability
    dimnames(cw)[[1]][1] <- "wrong"
    expect_error(setFishing(params, catchability = cw),
                 "The gear dimnames in the catchability array do not match the gear names.")
    cw <- catchability
    dimnames(cw)[[2]][1] <- "wrong"
    expect_error(setFishing(params, catchability = cw),
                 "The species dimnames in the catchability array do not match the species names.")

    selectivity <- getSelectivity(params)
    sw <- selectivity
    dimnames(sw)[[1]][1] <- "wrong"
    expect_error(setFishing(params, selectivity = sw),
                 "The gear dimnames in the selectivity array do not match the gear names.")
    sw <- selectivity
    dimnames(sw)[[2]][1] <- "wrong"
    expect_error(setFishing(params, selectivity = sw),
                 "The species dimnames in the selectivity array do not match the species names.")
    sw <- selectivity
    dimnames(sw)[[3]][1] <- "wrong"
    expect_warning(setFishing(params, selectivity = sw),
                 "I have changed the size dimnames in the selectivity array to agree with mizer conventions.")
})

test_that("Dimensions after number of gears has increased", {
    params <- NS_params_small
    gear_params(params)$gear <- params@species_params$species
    no_gears <- nrow(params@species_params)
    expect_identical(dim(params@catchability)[[1]], no_gears)
    expect_identical(dim(params@selectivity)[[1]], no_gears)
    expect_identical(dimnames(params@catchability)[[1]],
                     params@species_params$species)
    expect_identical(dimnames(params@selectivity)[[1]],
                     params@species_params$species)
    # The initial effort has also changed to the edition-specific default
    effort <- rep(ifelse(defaults_edition() < 2, 0, 1), no_gears)
    names(effort) <- params@species_params$species
    expect_identical(params@initial_effort, effort)
})

test_that("Duplicate gear-species pairs give error", {
    gp <- rbind(NS_params_small@gear_params, NS_params_small@gear_params[1, ])
    expect_error(gear_params(params) <- gp,
                 "Some species - gear pairs appear more than once")
})

test_that("Non-existing species give error", {
    gp <- NS_params_small@gear_params
    gp$species[[1]] <- "test"
    expect_error(gear_params(params) <- gp,
                 "The gear_params dataframe contains species that do not exist in the model.")
})

test_that("Can get and set selectivity slot", {
    params <- NS_params_small
    new <- 2 * selectivity(params)
    comment(new) <- "test"
    selectivity(params) <- new
    expect_identical(selectivity(params), new)
    expect_identical(getSelectivity(params), new)
})
test_that("Can get and set catchability slot", {
    params <- NS_params_small
    new <- 2 * catchability(params)
    comment(new) <- "test"
    catchability(params) <- new
    expect_identical(catchability(params), new)
    expect_identical(getCatchability(params), new)
})
test_that("Can get and set initial_effort slot", {
    params <- NS_params_small
    expect_identical(getInitialEffort(params), initial_effort(params))
    new <- 2 * initial_effort(params)
    initial_effort(params) <- new
    expect_identical(initial_effort(params), new)
    expect_identical(getInitialEffort(params), new)
})

# Bin-averaged (second-order) selectivity ----
test_that("Default path leaves selectivity point-sampled", {
    params <- NS_params_small
    # The default flag is off, so calc_selectivity must reproduce a direct
    # point evaluation of the selectivity functions at params@w.
    expect_false(second_order_w(params)[["bin_average"]])
    sel <- calc_selectivity(params)
    expect_identical(sel, params@selectivity[])
})

test_that("Bin-averaging changes selectivity and persists through setParams", {
    params <- NS_params_small
    sel_point <- calc_selectivity(params)
    second_order_w(params) <- c(bin_average = TRUE)
    # The flag is recorded
    expect_true(second_order_w(params)[["bin_average"]])
    sel_avg <- calc_selectivity(params)
    expect_equal(dim(sel_avg), dim(sel_point))
    # It really differs from the point-sampled version
    expect_gt(max(abs(sel_avg - sel_point)), 1e-3)
    # The setParams pipeline (re-run by the setter) stored the bin average
    expect_equal(params@selectivity[], sel_avg)
    # Still a valid selectivity (in [0, 1])
    expect_true(all(sel_avg >= 0 & sel_avg <= 1))
})

test_that("Knife-edge bin average is the analytic partial-bin fraction", {
    params <- NS_params_small
    w <- params@w
    dw <- params@dw
    j <- 8
    # Knife edge strictly inside bin j
    edge <- w[j] + 0.37 * dw[j]
    gp <- gear_params(params)
    gp$sel_func <- "knife_edge"
    gp$knife_edge_size <- edge
    gear_params(params) <- gp
    second_order_w(params) <- c(bin_average = TRUE)
    sel <- calc_selectivity(params)["Pelagic", "Herring", ]
    # Zero below the straddling bin, one above it
    expect_true(all(sel[seq_len(j - 1)] == 0))
    expect_true(all(abs(sel[(j + 1):length(sel)] - 1) < 1e-12))
    # In the straddling bin: fraction of the (linear-w) bin above the edge
    wjp1 <- w[j] + dw[j]
    frac <- (wjp1 - edge) / (wjp1 - w[j])
    # Composite-midpoint rule converges to the analytic fraction
    expect_equal(sel[j], frac, tolerance = 1e-2, ignore_attr = TRUE)
})

test_that("Smooth selectivity bin average matches a fine reference", {
    params <- NS_params_small
    gp <- gear_params(params)
    gp$sel_func <- "sigmoid_weight"
    gp$sigmoidal_weight <- 50
    gp$sigmoidal_sigma <- 3
    gear_params(params) <- gp
    second_order_w(params) <- c(bin_average = TRUE)
    sel <- calc_selectivity(params)["Pelagic", "Herring", ]
    w <- params@w
    dw <- params@dw
    ref <- vapply(seq_along(w), function(jj) {
        wjp1 <- w[jj] + dw[jj]
        ww <- seq(w[jj], wjp1, length.out = 4000)
        Sv <- sigmoid_weight(ww, sigmoidal_weight = 50, sigmoidal_sigma = 3)
        mean((Sv[-1] + Sv[-length(Sv)]) / 2)
    }, numeric(1))
    expect_equal(sel, ref, tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("Projection with bin-averaged selectivity runs finite", {
    params <- NS_params_small
    second_order_w(params) <- c(bin_average = TRUE)
    sim <- project(params, t_max = 2)
    expect_true(all(is.finite(sim@n)))
})
