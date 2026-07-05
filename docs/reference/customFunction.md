# Replace a mizer function with a custom version

**\[experimental\]** This function allows you to make arbitrary changes
to how mizer works by allowing you to replace any mizer function with
your own version. You should do this only as a last resort, when you
find that you can not use the standard mizer extension mechanism to
achieve your goal.

## Usage

``` r
customFunction(name, fun)
```

## Arguments

- name:

  Name of mizer function to replace

- fun:

  The custom function to use as replacement

## Value

No return value, called for side effects

## Details

If the function you need to overwrite is one of the mizer rate
functions, then you should use
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
instead of this function. Similarly you should use
`resource_dynamics()<-` to change the resource dynamics and
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md)
to change the density-dependence in reproduction. You should also
investigate whether you can achieve your goal by introducing additional
ecosystem components with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md).

If you find that your goal really does require you to overwrite a mizer
function, please also create an issue on the mizer issue tracker at
<https://github.com/sizespectrum/mizer/issues> to describe your goal,
because it will be interesting to the mizer community and may motivate
future improvements to the mizer functionality.

Note that `customFunction()` only overwrites the function used by the
mizer code. It does not overwrite the function that is exported by
mizer. This will become clear when you run the code in the Examples
section.

This function does not in any way check that your replacement function
is compatible with mizer. Calling this function can totally break mizer.
However you can always undo the effect by reloading mizer with

    detach(package:mizer, unload = TRUE)
    library(mizer)

## See also

"Extending mizer":
[`vignette("extending-mizer", package = "mizer")`](https://sizespectrum.org/mizer/articles/extending-mizer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fake_project <- function(...) "Fake"
customFunction("project", fake_project)
mizer::project(NS_params) # This will print "Fake"
project(NS_params) # This will still use the old project() function
# To undo the effect:
customFunction("project", project)
mizer::project(NS_params) # This will again use the old project()
} # }
```
