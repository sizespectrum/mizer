test_that("setMetadata works", {
    params <- NS_params
    params <- setMetadata(params, title = "title", 
                          description = "description",
                          authors = "Gustav Delius",
                          new = "new")
    # This should not change time_modified
    expect_identical(params@time_modified, NS_params@time_modified)
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
