# Density dependence used in these tests. The functions need to be found by
# name from inside mizer, so they live in the global environment.
assign("test_immigrationRDD", function(rdi, species_params, params, ...) {
    BevertonHoltRDD(rdi, species_params) + other_params(params)$immigration
}, envir = globalenv())
# Couples the species: all of them share one Beverton-Holt limit
assign("test_pooledRDD", function(rdi, species_params, ...) {
    rdi / (1 + sum(rdi) / species_params$R_max)
}, envir = globalenv())
assign("test_seasonalRDD", function(rdi, species_params, t, ...) {
    rdi * (1 + 0.1 * sin(2 * pi * t))
}, envir = globalenv())
assign("test_fixedRDD", function(rdi, species_params, ...) {
    rep(1e10, length(rdi))
}, envir = globalenv())

# The ratio of the reproduction rate to the rate required for steady state
rdd_ratio <- function(params) {
    unname(getRDD(params) / getRequiredRDD(params))
}

# NS_params_small produces far more eggs than it needs, so that its
# Beverton-Holt reproduction is saturated. Balanced at a lower reproduction
# level it leaves room for the reproduction to respond.
balanced_params <- setBevertonHolt(NS_params_small, reproduction_level = 0.25)

# A model with a fifth of its required reproduction coming from immigration,
# balanced by hand
with_immigration <- function(params = balanced_params) {
    other_params(params)$immigration <- 0.2 * unname(getRequiredRDD(params))
    target <- getRequiredRDD(params) - other_params(params)$immigration
    rdi <- getRDI(params)
    params@species_params$R_max <- unname(rdi * target / (rdi - target))
    params
}

# An extension that adds the immigration through a projectRDD() method
immigration_extension <- function(class, balanced) {
    registerS3method("projectRDD", class, function(params, rdi, ...) {
        NextMethod() + other_params(params)$immigration
    }, envir = asNamespace("mizer"))
    if (balanced) {
        registerS3method("balanceRDD", class, function(
            params, rdi = getRDI(params), rdd_target = getRequiredRDD(params),
            reproduction_level = NULL, ...) {
            immigration <- other_params(params)$immigration
            if (!is.null(reproduction_level)) {
                reproduction_level <- 1 - (1 - reproduction_level) *
                    rdd_target / (rdd_target - immigration)
            }
            NextMethod(rdd_target = rdd_target - immigration,
                       reproduction_level = reproduction_level)
        }, envir = asNamespace("mizer"))
    }
    params <- with_immigration()
    params@extensions <- setNames(NA_character_, class)
    coerceToExtensionClass(params)
}

perturbed <- function(params) {
    initialN(params) <- initialN(params) * 1.3
    params
}

# Beverton-Holt ----
test_that("balanceRDD() on Beverton-Holt is setBevertonHolt()", {
    params <- perturbed(NS_params_small)
    p1 <- balanceRDD(params)
    p2 <- setBevertonHolt(params)
    expect_equal(p1@species_params$R_max, p2@species_params$R_max)
    expect_equal(p1@species_params$erepro, params@species_params$erepro)
    expect_equal(rdd_ratio(p1), rep(1, 3))

    p1 <- suppressWarnings(balanceRDD(params, reproduction_level = 0.3))
    p2 <- suppressWarnings(setBevertonHolt(params, reproduction_level = 0.3))
    expect_equal(p1@species_params$R_max, p2@species_params$R_max)
    expect_equal(p1@species_params$erepro, p2@species_params$erepro)
})

test_that("balanceRDD() leaves alone species without a target", {
    params <- perturbed(NS_params_small)
    target <- getRequiredRDD(params)
    target[c(1, 3)] <- NA
    p <- balanceRDD(params, rdd_target = target)
    expect_identical(p@species_params$R_max[c(1, 3)],
                     params@species_params$R_max[c(1, 3)])
    expect_equal(rdd_ratio(p)[2], 1)
    expect_identical(balanceRDD(params, rdd_target = rep(NA, 3)), params)
})

test_that("balanceRDD() checks the length of its vectors", {
    expect_error(balanceRDD(NS_params_small, rdd_target = 1:2),
                 "one value for each species")
})

# Other built-in density dependence ----
test_that("balance_RickerRDD() keeps erepro and sets ricker_b", {
    params <- setReproduction(NS_params_small, RDD = "RickerRDD")
    p <- balanceRDD(params)
    expect_identical(p@species_params$erepro, params@species_params$erepro)
    expect_true(all(p@species_params$ricker_b >= 0))
    expect_equal(rdd_ratio(p), rep(1, 3))
    expect_equal(p@given_species_params$ricker_b, p@species_params$ricker_b)

    p <- balanceRDD(params, erepro = c(Cod = 0.5))
    expect_equal(p@species_params$erepro[3], 0.5, ignore_attr = TRUE)
    expect_identical(p@species_params$erepro[1:2],
                     params@species_params$erepro[1:2])
    expect_equal(rdd_ratio(p)[3], 1)
})

test_that("balance_RickerRDD() increases an erepro that is too small", {
    params <- setReproduction(NS_params_small, RDD = "RickerRDD")
    expect_warning(p <- balanceRDD(params, erepro = 1e-10),
                   "`erepro` has been increased to the smallest possible")
    expect_equal(p@species_params$ricker_b, rep(0, 3), ignore_attr = TRUE)
    expect_equal(rdd_ratio(p), rep(1, 3))
})

test_that("balance_SheperdRDD() keeps erepro and sheperd_c", {
    params <- setReproduction(NS_params_small, RDD = "SheperdRDD")
    params@species_params$sheperd_c <- 2
    p <- balanceRDD(params)
    expect_identical(p@species_params$erepro, params@species_params$erepro)
    expect_identical(p@species_params$sheperd_c, rep(2, 3),
                     ignore_attr = TRUE)
    expect_equal(rdd_ratio(p), rep(1, 3))

    params@species_params$sheperd_c <- NULL
    expect_error(balanceRDD(params), "must contain a column sheperd_c")
})

test_that("balance_noRDD() sets erepro", {
    params <- setRateFunction(perturbed(NS_params_small), "RDD", "noRDD")
    expect_no_message(p <- balanceRDD(params))
    expect_equal(rdd_ratio(p), rep(1, 3))
    expect_equal(p@given_species_params$erepro, p@species_params$erepro)
})

test_that("balance_constantRDD() sets constant_reproduction", {
    params <- NS_params_small
    params@rates_funcs$RDD <- "constantRDD"
    p <- balanceRDD(params)
    expect_equal(p@species_params$constant_reproduction,
                 getRequiredRDD(params), ignore_attr = TRUE)
    expect_equal(rdd_ratio(p), rep(1, 3))
})

test_that("a balancing function rejects arguments it does not accept", {
    params <- setReproduction(NS_params_small, RDD = "RickerRDD")
    expect_error(balanceRDD(params, reproduction_level = 0.5),
                 "`balance_RickerRDD\\(\\)` does not accept the argument `reproduction_level`")
    # A `NULL` argument has not been given
    expect_no_error(balanceRDD(params, reproduction_level = NULL))
})

# Custom density dependence ----
test_that("balanceRDD() solves for erepro when there is no balancing function", {
    params <- perturbed(with_immigration())
    params <- setRateFunction(params, "RDD", "test_immigrationRDD")
    expect_message(p <- balanceRDD(params),
                   "no `balance_test_immigrationRDD\\(\\)` function.*adjusting `erepro` numerically")
    expect_equal(rdd_ratio(p), rep(1, 3))
    expect_identical(p@species_params$R_max, params@species_params$R_max)
    expect_equal(p@given_species_params$erepro, p@species_params$erepro)
    expect_no_message(balanceRDD(params, info_level = 0))
})

test_that("balanceRDD() solves numerically for density dependence that couples species", {
    params <- setRateFunction(perturbed(NS_params_small), "RDD",
                              "test_pooledRDD")
    params@species_params$R_max <- 10 * sum(getRDI(params))
    p <- balanceRDD(params, info_level = 0)
    expect_equal(rdd_ratio(p), rep(1, 3))
})

test_that("only erepro can be adjusted without a balancing function", {
    params <- setRateFunction(with_immigration(), "RDD", "test_immigrationRDD")
    expect_error(balanceRDD(params, reproduction_level = 0.5),
                 "can only balance the reproduction by adjusting `erepro`")
})

test_that("balanceRDD() warns when no erepro gives the required rate", {
    params <- setRateFunction(NS_params_small, "RDD", "test_fixedRDD")
    expect_warning(p <- balanceRDD(params),
                   "no value of `erepro` gives the required reproduction rate")
    expect_identical(p@species_params, params@species_params)
})

test_that("balanceRDD() refuses density dependence that depends on time", {
    params <- setRateFunction(NS_params_small, "RDD", "test_seasonalRDD")
    expect_warning(p <- balanceRDD(params), "depends on time")
    expect_identical(p@species_params, params@species_params)
})

test_that("balanceRDD() stops for a present species without egg production", {
    params <- NS_params_small
    other_params(params)$immigration <- rep(0, 3)
    params <- setRateFunction(params, "RDD", "test_immigrationRDD")
    params@species_params$erepro[1] <- 0
    expect_error(balanceRDD(params, info_level = 0),
                 "no reproduction.*: Sprat\\.")
})

# Extensions ----
test_that("balanceRDD() balances through a projectRDD() method without a balanceRDD() method", {
    params <- perturbed(immigration_extension("mizerTestRDDUnbalanced",
                                              balanced = FALSE))
    expect_identical(params@rates_funcs$RDD, "BevertonHoltRDD")
    expect_message(p <- balanceRDD(params),
                   "`mizerTestRDDUnbalanced` changes the reproduction rate")
    expect_equal(rdd_ratio(p), rep(1, 3))
    # Beverton-Holt parameters can not be preserved through it
    expect_error(balanceRDD(params, reproduction_level = 0.5),
                 "can only balance the reproduction by adjusting `erepro`")
})

test_that("a balanceRDD() method passes the target on with NextMethod()", {
    params <- perturbed(immigration_extension("mizerTestRDDBalanced",
                                              balanced = TRUE))
    expect_no_message(p <- balanceRDD(params))
    expect_equal(rdd_ratio(p), rep(1, 3))
    expect_identical(p@species_params$erepro, params@species_params$erepro)

    reproduction_level(p) <- 0.6
    expect_equal(reproduction_level(p), rep(0.6, 3), ignore_attr = TRUE,
                 tolerance = 1e-6)
    expect_equal(rdd_ratio(p), rep(1, 3))
})

test_that("balanceRDD() gives up when only some extensions take part", {
    immigration_extension("mizerTestRDDOuter", balanced = TRUE)
    immigration_extension("mizerTestRDDInner", balanced = FALSE)
    params <- with_immigration()
    params@extensions <- setNames(rep(NA_character_, 2),
                                  c("mizerTestRDDOuter", "mizerTestRDDInner"))
    params <- coerceToExtensionClass(params)
    expect_warning(p <- balanceRDD(params),
                   "while other extensions do have one")
    expect_identical(p@species_params, params@species_params)
})

# Reproduction level ----
test_that("rdd_elasticity_level() agrees with the Beverton-Holt level", {
    expect_equal(rdd_elasticity_level(NS_params_small),
                 reproduction_level(NS_params_small), tolerance = 1e-8)
})

test_that("the second-order scheme gives a balanced model too", {
    params <- NS_params_small
    second_order_w(params) <- TRUE
    params <- setReproduction(perturbed(params), RDD = "RickerRDD")
    p <- balanceRDD(params)
    expect_equal(rdd_ratio(p), rep(1, 3))
})
