test_that("setMetadata works", {
    params <- NS_params_small
    params <- setMetadata(params, title = "title", 
                          description = "description",
                          authors = "Gustav Delius",
                          new = "new")
    # This should not change time_modified
    expect_identical(params@time_modified, NS_params_small@time_modified)
    metadata <- getMetadata(params)
    expect_identical(metadata$title, "title")
    expect_identical(metadata$description, "description")
    expect_identical(metadata$authors, "Gustav Delius")
    expect_identical(metadata$new, "new")
    
    params <- setMetadata(params, title = "new",
                          new = NULL)
    metadata <- getMetadata(params)
    expect_identical(metadata$title, "new")
    expect_null(metadata$new)
})

test_that("getMetadata always includes automatic fields", {
    metadata <- getMetadata(NS_params_small)
    expect_identical(metadata$mizer_version, NS_params_small@mizer_version)
    expect_identical(metadata$extensions, NS_params_small@extensions)
    expect_identical(metadata$time_created, NS_params_small@time_created)
    expect_identical(metadata$time_modified, NS_params_small@time_modified)
})

test_that("setMetadata ignores automatic fields supplied in dots", {
    expect_message(
        params <- setMetadata(NS_params_small, mizer_version = "bad", time_created = 0),
        "set automatically by mizer"
    )
    metadata <- getMetadata(params)
    expect_identical(metadata$mizer_version, NS_params_small@mizer_version)
    expect_identical(metadata$time_created, NS_params_small@time_created)
})
