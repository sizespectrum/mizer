species_params <- NS_species_params
w <- NS_params@w

test_that("sigmoid_length works", {
    expect_error(sigmoid_length(w, 20, 30, species_params = species_params),
                 "The selectivity function needs the weight-length parameters ")
    species_params$a <- 0.5
    species_params$b <- 3
    expect_length(sigmoid_length(w, 20, 30, species_params = species_params[1, ]),
                  length(w))
})

test_that("double_sigmoid_length works", {
    expect_error(double_sigmoid_length(w, 20, 30, 40, 50, 
                                       species_params = species_params),
                 "The selectivity function needs the weight-length parameters ")
    species_params$a <- 0.5
    species_params$b <- 3
    expect_length(double_sigmoid_length(w, 20, 30, 40, 50,
                                        species_params = species_params[1, ]),
                  length(w))
    expect_error(double_sigmoid_length(w, 20, 30, 40, 30,
                                        species_params = species_params[1, ]),
                 "l50_right not less than l25_right")
})

test_that("knife_edge works", {
    expect_length(knife_edge(w, 20, species_params = species_params[1, ]),
                  length(w))
})

test_that("sigmoid_weight works", {
    expect_length(knife_edge(w, 20, 2), length(w))
})
