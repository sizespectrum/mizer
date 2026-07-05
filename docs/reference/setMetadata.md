# Set metadata for a model

**\[experimental\]** Setting metadata is particularly important for
sharing your model with others. All metadata fields are optional and you
can also add other fields of your own choosing. If you set a value for a
field that already existed, the old value will be overwritten.

## Usage

``` r
setMetadata(
  params,
  title = NULL,
  description = NULL,
  authors = NULL,
  url = NULL,
  doi = NULL,
  ...
)

getMetadata(params)
```

## Arguments

- params:

  The MizerParams object for the model

- title:

  A string with the title for the model

- description:

  A string with a description of the model. This could for example
  contain information about any publications using the model.

- authors:

  An author entry or a list of author entries, where each author entry
  could either be just a name or could itself be a list with fields like
  `name`, `orcid`, possibly `email`.

- url:

  A URL where more information about the model can be found. This could
  be a blog post on the mizer blog, for example.

- doi:

  The digital object identifier for your model. To create a doi you can
  use online services like https://zenodo.org/ or https://figshare.com.

- ...:

  Additional metadata fields that you would like to add

## Value

`setMetadata()`: The MizerParams object with updated metadata

`getMetadata()`: A list with all metadata entries that have been set,
including at least `mizer_version`, `extensions`, `time_created` and
`time_modified`.

## Details

In addition to the metadata fields you can set by hand, there are four
fields that are set automatically by mizer:

- `mizer_version` The version string of the mizer version under which
  the model was created or last upgraded. Can be compared to the current
  version which is obtained with `packageVersion("mizer")`. The purpose
  of this field is that if the model is not working as expected in the
  current version of mizer, you can go back to the older version under
  which presumably it was working.

- `extensions` A named vector of strings where each name is the name of
  and extension package needed to run the model and each value is a
  string giving the information that the remotes package needs to
  install the correct version of the extension package, see
  https://remotes.r-lib.org/. This field is set by the extension
  packages.

- `time_created` A POSIXct date-time object with the creation time.

- `time_modified` A POSIXct date-time object with the last modified
  time.

Setting the metadata with this function does not count as a modification
of the object, so the `time_modified` field will not be updated.

## Examples

``` r
params <- setMetadata(NS_params,
    title = "North Sea model",
    description = "A multi-species model of the North Sea fish community.",
    authors = list(list(name = "Finlay Scott", email = "finlay@example.com")))
getMetadata(params)$title
#> [1] "North Sea model"
getMetadata(params)$authors[[1]]$name
#> [1] "Finlay Scott"
```
