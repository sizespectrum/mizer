test_that(".onAttach shows a message the first time a version is seen", {
    config_dir <- withr::local_tempdir()
    withr::local_envvar(R_USER_CONFIG_DIR = config_dir)
    current <- as.character(utils::packageVersion("mizer"))
    version_file <- file.path(config_dir, "R", "mizer", "last_seen_version")

    expect_message(mizer:::.onAttach(), current)
    expect_true(file.exists(version_file))
    expect_equal(readLines(version_file), current)
})

test_that(".onAttach stays silent on a repeat load of the same version", {
    config_dir <- withr::local_tempdir()
    withr::local_envvar(R_USER_CONFIG_DIR = config_dir)

    suppressMessages(mizer:::.onAttach())
    expect_no_message(mizer:::.onAttach())
})

test_that(".onAttach shows the message again after a version change", {
    config_dir <- withr::local_tempdir()
    withr::local_envvar(R_USER_CONFIG_DIR = config_dir)
    version_dir <- file.path(config_dir, "R", "mizer")
    dir.create(version_dir, recursive = TRUE)
    writeLines("0.0.1", file.path(version_dir, "last_seen_version"))

    current <- as.character(utils::packageVersion("mizer"))
    expect_message(mizer:::.onAttach(), current)
})

test_that(".onAttach never errors even if the config path can't be written", {
    blocker <- withr::local_tempfile()
    writeLines("not a directory", blocker)
    withr::local_envvar(R_USER_CONFIG_DIR = blocker)

    expect_no_error(suppressMessages(mizer:::.onAttach()))
})
