# Often I need to test that a MizerParams or MizerSim object has not
# changed except for the time_modified
expect_unchanged <- function(object, expected) {
    if (is(object, "MizerParams")) {
        object@time_modified <- expected@time_modified
    }
    if (is(object, "MizerSim")) {
        object@params@time_modified <- expected@params@time_modified
    }
    
    expect_equal(object, expected)
}