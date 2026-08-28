# Get the environment holding mizer's dynamic marker classes

The environment is attached to the search path, lazily, the first time a
dispatching extension needs a dynamic marker class. A session using base
mizer, or only metadata-only extensions, therefore never acquires the
extra search path entry.

## Usage

``` r
extensionClassEnvironment(create = TRUE)
```

## Arguments

- create:

  Logical. If `TRUE`, create and attach the environment when it does not
  yet exist. If `FALSE`, return `NULL` instead.

## Value

The attached environment, or `NULL` when it does not exist and
`create = FALSE`.

## Details

The search path is where this metadata has to live. `.GlobalEnv` does
not work, because `R CMD check` empties it with `cleanEx()` between
examples and users clear it themselves (#587). An environment that mizer
merely holds privately does not work either, however well it survives: a
params object saved before mizer started creating marker classes
dynamically carries a class attribute whose package slot names the
extension package, and `methods` then resolves the class with an
inheriting lookup that starts in that package's namespace and runs on
through `.GlobalEnv` and the rest of the search path. It never passes
through mizer's own namespace, so only the search path catches every
case.

`cleanEx()` also detaches whatever an example added to the search path,
so the environment lasts a whole check run only when it was attached
while the extension package was being loaded, before `R CMD check`
recorded the search path. That is what happens when an extension
registers itself from `.onLoad`. Losing it the other way is not fatal:
[`markerClassPresent()`](https://sizespectrum.org/mizer/reference/markerClassPresent.md)
then reports the classes as gone and the next registration rebuilds
them.

The environment is deliberately left mutable, because rebuilding an
extension chain requires removing and re-parenting its classes. The S4
package identity of the classes stays `.GlobalEnv`, which is what marks
them as mutable dynamic classes rather than static classes owned by an
extension namespace.
