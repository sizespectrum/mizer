context("Summary methods")

test_that("get_size_range_array",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
    size_n <- get_size_range_array(params)
    expect_that(all(size_n), is_true())
    size_n <- get_size_range_array(params, min_w = 1)
    expect_that(!all(size_n[,which(params@w < 1)]), is_true())
    expect_that(all(size_n[,which(params@w >= 1)]), is_true())
    size_n <- get_size_range_array(params, max_w = 100)
    expect_that(all(size_n[,which(params@w <= 100)]), is_true())
    expect_that(!all(size_n[,which(params@w > 1)]), is_true())
    size_n <- get_size_range_array(params, min_w = 1, max_w = 100)
    expect_that(!all(size_n[,which(params@w > 100)]), is_true())
    expect_that(!all(size_n[,which(params@w < 1)]), is_true())
    expect_that(all(size_n[,which((params@w >= 1) & (params@w<=100))]), is_true())
    size_n <- get_size_range_array(params, min_l = 1)
    min_w <- params@species_params$a * 1 ^ params@species_params$b
    for (sp in 1:nrow(params@species_params)){ 
	expect_that(all(size_n[sp,which(params@w >= min_w[sp])]), is_true())
	expect_that(!all(size_n[sp,which(params@w < min_w[sp])]), is_true())
    }
    size_n <- get_size_range_array(params, max_l = 100)
    max_w <- params@species_params$a * 100 ^ params@species_params$b
    for (sp in 1:nrow(params@species_params)){ 
	expect_that(all(size_n[sp,which(params@w <= max_w[sp])]), is_true())
	expect_that(!all(size_n[sp,which(params@w > max_w[sp])]), is_true())
    }
    size_n <- get_size_range_array(params, min_l = 1, max_l = 100)
    min_w <- params@species_params$a * 1 ^ params@species_params$b
    max_w <- params@species_params$a * 100 ^ params@species_params$b
    for (sp in 1:nrow(params@species_params)){ 
	expect_that(all(size_n[sp,which((params@w <= max_w[sp]) & (params@w >= min_w[sp]))]), is_true())
	expect_that(!all(size_n[sp,which(params@w < min_w[sp])]), is_true())
	expect_that(!all(size_n[sp,which(params@w > max_w[sp])]), is_true())
    }
    size_n <- get_size_range_array(params, min_w = 1, max_l = 100)
    min_w <- rep(1,nrow(params@species_params))
    max_w <- params@species_params$a * 100 ^ params@species_params$b
    for (sp in 1:nrow(params@species_params)){ 
	expect_that(all(size_n[sp,which((params@w <= max_w[sp]) & (params@w >= min_w[sp]))]), is_true())
	expect_that(!all(size_n[sp,which(params@w < min_w[sp])]), is_true())
	expect_that(!all(size_n[sp,which(params@w > max_w[sp])]), is_true())
    }
    size_n <- get_size_range_array(params, min_l = 1, max_w = 100)
    min_w <- params@species_params$a * 1 ^ params@species_params$b
    max_w <- rep(100,nrow(params@species_params))
    for (sp in 1:nrow(params@species_params)){ 
	expect_that(all(size_n[sp,which((params@w <= max_w[sp]) & (params@w >= min_w[sp]))]), is_true())
	expect_that(!all(size_n[sp,which(params@w < min_w[sp])]), is_true())
	expect_that(!all(size_n[sp,which(params@w > max_w[sp])]), is_true())
    }
})



