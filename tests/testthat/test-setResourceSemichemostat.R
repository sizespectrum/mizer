test_that("setResourceSemichemostat works", {
    # semichemostat rate
    sc <- function(params) {
        r <- getResourceRate(params)
        c <- getResourceCapacity(params)
        N <- initialNResource(params)
        mu <- getResourceMort(params)
        N_steady <- r * c / (mu + r)
        sel <- (mu + r) > 0
        all.equal(N_steady[sel], N[sel])
    }
    
    params <- NS_params
    
    # setting rate
    rate <- getResourceMort(params)
    params <- setResourceSemichemostat(params, resource_rate = rate)
    expect_identical(getResourceRate(params), rate)
    expect_true(sc(params))
    # single value
    p1 <- setResourceSemichemostat(params, resource_rate = 2)
    expect_equal(unname(getResourceRate(p1)), 
                     2 * params@w_full ^ (p1@resource_params$n - 1))
    
    # setting capacity
    capacity <- getResourceCapacity(params) * 2
    p1 <- setResourceSemichemostat(params, resource_capacity = capacity)
    expect_identical(getResourceCapacity(p1), capacity)
    expect_true(sc(p1))
    # single value
    p1 <- setResourceSemichemostat(params, resource_capacity = 3e11)
    expect_equal(getResourceCapacity(p1)[1],
                     3e11 * params@w_full[1] ^ (-p1@resource_params$lambda))
    
    # setting level
    level <- getResourceLevel(params) / 5
    p1 <- setResourceSemichemostat(params, resource_level = level)
    expect_equal(getResourceLevel(p1), level)
    expect_true(sc(p1))
    # single value
    p1 <- setResourceSemichemostat(params, resource_level = 1/3)
    expect_equal(getResourceCapacity(p1), initialNResource(p1) * 3)
    # illegal value
    expect_error(setResourceSemichemostat(params, resource_level = 1),
                 "The 'resource_level' must always be")
    
    # Without resource argument
    params@initial_n <- params@initial_n * 20
    p1 <- setResourceSemichemostat(params)
    expect_true(sc(p1))
})
