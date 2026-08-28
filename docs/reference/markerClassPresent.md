# Is a marker class definition backed by a live binding?

[`methods::isClass()`](https://rdrr.io/r/methods/findClass.html) can
keep returning `TRUE` after [`rm()`](https://rdrr.io/r/base/rm.html) has
deleted the class metadata binding. Package-owned classes are not
susceptible to that global-environment wipe. For a dynamically defined
class, however, the binding must still exist either in mizer's
persistent class environment or, for compatibility with classes created
by earlier mizer versions, in `.GlobalEnv`.

## Usage

``` r
markerClassPresent(class)
```

## Arguments

- class:

  Character string — the S4 class name to test.

## Value

`TRUE` if the class can still be resolved.
