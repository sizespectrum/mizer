test_that("setResourceLogistic works", {
    # logistic rate
    sc <- function(params) {
        r <- getResourceRate(params)
        c <- getResourceCapacity(params)
        N <- initialNResource(params)
        mu <- getResourceMort(params)
        N_steady <- c * (r - mu) / r
        sel <- (r - mu) > 0
        all.equal(N_steady[sel], N[sel])
    }
    
    params <- NS_params
    
    # setting rate
    rate <- getResourceMort(params) * 1.1
    params <- setResourceLogistic(params, resource_rate = rate)
    expect_identical(getResourceRate(params), rate)
    expect_true(sc(params))
    # single value
    p1 <- setResourceLogistic(params, resource_rate = 30)
    expect_equal(unname(getResourceRate(p1)), 
                     30 * params@w_full ^ (p1@resource_params$n - 1))
    
    # setting capacity
    capacity <- getResourceCapacity(params) * 2
    p1 <- setResourceLogistic(params, resource_capacity = capacity)
    expect_identical(getResourceCapacity(p1), capacity)
    expect_true(sc(p1))
    # single value
    p1 <- setResourceLogistic(params, resource_capacity = 3e11)
    expect_equal(getResourceCapacity(p1)[1],
                     3e11 * params@w_full[1] ^ (-p1@resource_params$lambda))
    # make sure cutoff is implemented
    sel <- params@w_full > params@resource_params$w_pp_cutoff
    expect_equal(unname(getResourceCapacity(p1)[sel]),
                 rep(0, sum(sel)))
    
    # setting level
    level <- getResourceLevel(params) / 5
    p1 <- setResourceLogistic(params, resource_level = level)
    expect_equal(getResourceLevel(p1), level)
    expect_true(sc(p1))
    # single value
    p1 <- setResourceLogistic(params, resource_level = 1/3)
    expect_equal(getResourceCapacity(p1), initialNResource(p1) * 3)
    # illegal value
    expect_error(setResourceLogistic(params, resource_level = 1),
                 "The 'resource_level' must always be")
    
    # Without resource argument
    params@initial_n <- params@initial_n / 20
    p1 <- setResourceLogistic(params)
    expect_true(sc(p1))
})
