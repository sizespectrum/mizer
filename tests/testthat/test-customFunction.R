test_that("'customFunction()' works", {
    sim <- project(NS_params, t_max = 1)
    project_original <- project
    fake_project <- function(...) "Fake"
    customFunction("project", fun = fake_project)
    expect_identical(mizer::project(NS_params), "Fake")
    # To undo the effect:
    customFunction("project", project_original)
    expect_identical(mizer::project(NS_params, t_max = 1), sim)
    # fun argument should be a function
    expect_error(customFunction("project", fun = 1),
                 "fun is not a function")
    expect_error(customFunction("projet", fun = fake_project),
                 "There is no mizer function called projet")
})
