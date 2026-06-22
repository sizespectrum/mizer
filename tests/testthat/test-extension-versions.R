test_that("upgrade is dispatched as an S3 generic with our methods", {
    methods <- as.character(utils::.S3methods("upgrade"))
    expect_true("upgrade.MizerParams" %in% methods)
    expect_true("upgrade.MizerSim" %in% methods)
    # The base R method must still be reachable.
    expect_false(is.null(getS3method("upgrade", "packageStatus", optional = TRUE)))
})

test_that("extensionRequirements / extensionVersions read both slot forms", {
    # legacy character form: requirements only, versions all NA
    char <- c(extA = "owner/repo", extB = NA_character_)
    expect_identical(extensionRequirements(char), char)
    expect_identical(extensionVersions(char),
                     setNames(c(NA_character_, NA_character_), c("extA", "extB")))

    # versioned list form
    lst <- makeExtensions(c(extA = "owner/repo"), c(extA = "1.2.3"))
    expect_identical(extensionRequirements(lst), c(extA = "owner/repo"))
    expect_identical(extensionVersions(lst), c(extA = "1.2.3"))

    # empty
    expect_identical(extensionRequirements(character()), character())
})

test_that("recordExtension stamps without disturbing other entries", {
    p <- NS_params
    p@extensions <- list(other = c(requirement = NA_character_,
                                   version = NA_character_),
                         extA = c(requirement = "owner/repo",
                                  version = NA_character_))
    before <- p@extensions

    # version = NULL on an existing entry is a no-op on the whole slot
    p2 <- recordExtension(p, "extA")
    expect_identical(p2@extensions, before)

    # stamping updates only the named entry
    p3 <- recordExtension(p, "extA", version = "2.0.0")
    expect_identical(p3@extensions$other, before$other)
    expect_identical(unname(p3@extensions$extA[["version"]]), "2.0.0")
})

test_that("recordExtension can stamp onto an empty slot", {
    p <- NS_params
    expect_length(p@extensions, 0)
    p <- recordExtension(p, "extA", version = "1.0.0")
    expect_true(is.list(p@extensions))
    expect_identical(unname(p@extensions$extA[["requirement"]]), NA_character_)
    expect_identical(unname(p@extensions$extA[["version"]]), "1.0.0")
})

test_that("registerExtensions accepts the versioned list form", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())
    ext <- paste0("mizerTestVer", Sys.getpid())
    lst <- makeExtensions(setNames(NA_character_, ext), character())
    expect_invisible(registerExtensions(lst))
    expect_identical(getRegisteredExtensions(), setNames(NA_character_, ext))
})

test_that("extension_needs_upgrading reacts to missing/stale stamps", {
    p <- NS_params
    # No extensions -> never needs upgrading on extension grounds
    expect_false(extension_needs_upgrading(p))

    # An extension whose package is not installed is ignored
    p@extensions <- makeExtensions(c(definitelyMissingPkg = "owner/repo"),
                                   character())
    expect_false(extension_needs_upgrading(p))

    # mizer itself, with a stale stamp, does need upgrading; with the installed
    # version it does not. (mizer is, of course, installed.)
    p@extensions <- makeExtensions(c(mizer = "x"), c(mizer = "0.0.1"))
    expect_true(extension_needs_upgrading(p))
    p@extensions <- makeExtensions(
        c(mizer = "x"),
        c(mizer = as.character(utils::packageVersion("mizer"))))
    expect_false(extension_needs_upgrading(p))
})
