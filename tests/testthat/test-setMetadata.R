test_that("setMetadata works", {
    params <- setMetadata(NS_params, title = "title", 
                          description = "description",
                          authors = "Gustav Delius",
                          orcid = c("Gustav Delius" = "0000-0003-4092-8228"))
    metadata <- getMetadata(params)
    expect_identical(metadata$title, "title")
    expect_identical(metadata$description, "description")
    expect_identical(metadata$authors, "Gustav Delius")
    expect_identical(metadata$orcid, c("Gustav Delius" = "0000-0003-4092-8228"))
    
    # multiple authors
    params <- setMetadata(params, authors = c("Gustav Delius", "Donald Duck"))
    expect_identical(getMetadata(params)$authors,
                     c("Gustav Delius", "Donald Duck"))
})
