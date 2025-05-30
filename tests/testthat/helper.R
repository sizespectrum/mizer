# Test that a MizerParams or MizerSim object has not changed except for the
# time_modified and perhaps a reordering of the species_params columns.
expect_unchanged <- function(object, expected) {
    if (is(object, "MizerParams")) {
        # has updated time_modified
        expect_false(identical(object@time_modified, expected@time_modified))
        object@time_modified <- expected@time_modified
        sp <- object@species_params
        sp_expected <- expected@species_params
    }
    if (is(object, "MizerSim")) {
        # has updated time_modified
        expect_false(identical(object@params@time_modified,
                               expected@params@time_modified))
        object@params@time_modified <- expected@params@time_modified
        sp <- object@params@species_params
        sp_expected <- expected@params@species_params
    }
    # Check that the species_params are unchanged except for a
    # reordering of the dataframe columns
    expect_equal(sp[, sort(names(sp))],
                     sp_expected[, sort(names(sp_expected))])
    # And if they were the same, remove the reordering
    if (is(object, "MizerParams")) {
        object@species_params <- sp_expected
    } else if (is(object, "MizerSim")) {
        object@params@species_params <- sp_expected
    }

    expect_equal(object, expected)
}
