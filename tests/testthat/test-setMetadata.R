test_that("setMetadata works", {
    params <- NS_params
    params <- setMetadata(params, title = "title", 
                          description = "description",
                          authors = "Gustav Delius",
                          new = "new")
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
    expect_identical(metadata$mizer_version, packageVersion("mizer"))
})
