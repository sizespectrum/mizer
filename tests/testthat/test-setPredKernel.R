params <- NS_params
no_sp <- nrow(params@species_params)

## setPredKernel ----
test_that("setPredKernel works", {
    expect_equal(setPredKernel(params), params)
    expect_equal(setPredKernel(params, pred_kernel = NULL), 
                     params)
    params@species_params$pred_kernel_type <- "box"
    params@species_params$ppmr_min <- 2
    expect_error(setPredKernel(params), 
                 "missing from the parameter dataframe: ppmr_max")
    params@species_params$ppmr_max <- 4
    p2 <- setPredKernel(params)
    pred_kernel <- 2 * getPredKernel(params)
    expect_error(setPredKernel(params, pred_kernel[1:2, ]),
                 "incorrect number of dimensions")
    expect_error(setPredKernel(params, pred_kernel - 1),
                 "pred_kernel >= 0 are not true")
    p2 <- setPredKernel(params, pred_kernel)
    expect_equivalent(p2@pred_kernel, pred_kernel)
    expect_identical(p2@pred_kernel, getPredKernel(p2))
})

test_that("Comment works on pred_kernel", {
    params <- NS_params
    # if no comment, it is set automatically
    pred_kernel <- getPredKernel(params)
    params <- setPredKernel(params, pred_kernel = pred_kernel)
    expect_identical(comment(params@pred_kernel), "set manually")
    
    # comment is stored
    comment(pred_kernel) <- "test"
    params <- setPredKernel(params, pred_kernel = pred_kernel)
    expect_identical(comment(params@pred_kernel), "test")
    
    # if no comment, previous comment is kept
    comment(pred_kernel) <- NULL
    params <- setPredKernel(params, pred_kernel = pred_kernel)
    expect_identical(comment(params@pred_kernel), "test")
    
    # no message when nothing changes
    expect_message(setPredKernel(params), NA)
    # but message when a change is not stored due to comment
    beta <- params@species_params$beta
    params@species_params$beta <- 1
    expect_message(setPredKernel(params),
                   "You have set a custom predation kernel")
    # Can reset
    params@species_params$beta <- beta
    p <- setPredKernel(params, reset = TRUE)
    expect_equal(p@pred_kernel, pred_kernel)
    expect_warning(setPredKernel(params, pred_kernel = pred_kernel,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

# getPredKernel ----
test_that("getPredKernel has correct dimnames",{
    pred_kernel <- getPredKernel(params)
    expect_identical(dimnames(pred_kernel)$sp, 
                     dimnames(params@initial_n)$sp)
    expect_identical(dimnames(pred_kernel)$w_pred, 
                     dimnames(params@initial_n)$w)
    expect_identical(dimnames(pred_kernel)$w_prey, 
                     as.character(signif(params@w_full, 3)))
})
test_that("getting and setting pred kernel leads to same dynamics" ,{
    params <- NS_params
    params <- setPredKernel(params, pred_kernel = getPredKernel(params))
    sim1 <- project(NS_params, t_max = 0.1)
    sim2 <- project(params, t_max = 0.1)
    expect_equal(finalN(sim1), finalN(sim2), tolerance = 1e-4)
})

test_that("Can get and set pred_kernel slot", {
    params <- NS_params
    new <- 2 * pred_kernel(params)
    comment(new) <- "test"
    pred_kernel(params) <- new
    expect_identical(pred_kernel(params), new)
})

## get_phi ----
test_that("get_phi works", {
    NS_species_params$pred_kernel_type <- "box"
    NS_species_params$ppmr_min <- 2
    NS_species_params$ppmr_max <- 4
    phi <- get_phi(NS_species_params, 1:5)
    expect_identical(phi[1, ], phi[2, ])
    expect_identical(phi[1, 1], 0)
    expect_identical(phi[1, 2], 1)
    expect_identical(phi[1, 5], 0)
})