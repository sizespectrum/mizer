test_that("Time resampling in project behaves as documented", {
    devtools::load_all(".")
    params <- NS_params
    gear_names <- unique(gear_params(params)$gear)

    # Create effort array with irregular times
    times <- c(1, 10, 11)
    effort <- array(1,
        dim = c(length(times), length(gear_names)),
        dimnames = list(time = times, gear = gear_names)
    )

    # Case 1: No t_max or t_save provided
    # Should preserve original times
    sim <- project(params, effort = effort, dt = 0.1)
    expect_equal(as.numeric(dimnames(sim@n)[[1]]), times)

    # Case 2: t_max provided
    # Should resample based on default t_save (which is 1 or inferred)
    # Here inferred t_save from first interval is 10 - 1 = 9.
    # So times would be 1, 10, 19...
    # Wait, the code says:
    # if (length(effort_times) > 1) {
    #     save_freq <- effort_times[2] - effort_times[1]
    # }
    # So save_freq becomes 9.
    # New times: seq(1, 20, by = 9) -> 1, 10, 19.
    # Time 11 is lost.
    sim <- project(params, effort = effort, t_max = 20, dt = 0.1)
    saved_times <- as.numeric(dimnames(sim@n)[[1]])

    expect_equal(saved_times, c(1, 10, 19))
    expect_false(11 %in% saved_times)

    # Case 3: t_save provided explicitly
    # t_save = 1. Should capture 1, 2, ..., 11, ...
    sim <- project(params, effort = effort, t_max = 20, t_save = 1, dt = 0.1)
    saved_times <- as.numeric(dimnames(sim@n)[[1]])
    expect_true(11 %in% saved_times)

    # Case 4: t_save provided but does not align with irregular time
    # t_save = 5. Times: 1, 6, 11, 16, 21...
    # Here 11 IS captured by coincidence (1 + 2*5 = 11).
    # Let's try t_save = 4. Times: 1, 5, 9, 13... 11 is lost.
    sim <- project(params, effort = effort, t_max = 20, t_save = 4, dt = 0.1)
    saved_times <- as.numeric(dimnames(sim@n)[[1]])
    expect_false(11 %in% saved_times)
})
