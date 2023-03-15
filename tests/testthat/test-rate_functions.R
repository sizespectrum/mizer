test_that("get_time_elements", {
    params <- NS_params
    sim <- project(params, effort = 1, t_max = 10, dt = 0.5, t_save = 0.5)
    expect_identical(get_time_elements(sim, as.character(3:4)),
                     get_time_elements(sim, 3:4))
    expect_identical(length(get_time_elements(sim, 3:4)),
                     dim(sim@n)[1])
    expect_identical(sum(get_time_elements(sim, 3:4)), 3L)
    expect_error(get_time_elements(sim, 3:50), 
                 "Time range is outside the time range of the model")
    expect_equivalent(which(get_time_elements(sim, seq(3, 4, by = 0.1))), 
                      c(7, 8, 9))
    # What if real years are used
    effort <- array(1, dim = c(19, 4),
                    dimnames = list(year = seq(1960, 1969, by = 0.5), 
                                    gear = c("Industrial", "Pelagic",
                                             "Otter", "Beam")
                                    )
                    )
    sim <- project(params, effort = effort, t_save = 0.5)
    expect_equivalent(which(get_time_elements(sim, 1965)), 11)
    expect_equivalent(which(get_time_elements(sim, "1965")), 11)
    expect_equivalent(which(get_time_elements(sim, 1965:1969)), 11:19)
})
