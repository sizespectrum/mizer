test_that("pred kernel functions", {
    ppmr <- NS_params@w_full / NS_params@w_full[1]
    phi <- lognormal_pred_kernel(ppmr, beta = 100, sigma = 2)
    expect_true(all(phi > 0))
    
    phit <- truncated_lognormal_pred_kernel(ppmr, beta = 100, sigma = 2)
    expect_equal(phi[2], phit[2])
    expect_equal(phit[61], 0)
    
    phib <- box_pred_kernel(ppmr = 1:5, ppmr_min = 2, ppmr_max = 4)
    expect_identical(phib[1], 0)
    expect_identical(phib[2], 1)
    expect_identical(phib[5], 0)
    
    phip <- power_law_pred_kernel(1:100, kernel_exp = 0,
                                  kernel_l_l = log(25) , kernel_u_l = 1000,
                                  kernel_l_r = log(75) , kernel_u_r = 1000)
    expect_equal(phip[10], 0)
    expect_equal(phip[30], 1)
    expect_equal(phip[90], 0)
})
