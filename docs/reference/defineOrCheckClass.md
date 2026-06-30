# Define an S4 class or verify it extends the expected parent

If `class` does not yet exist, defines it as a virtual-free S4 class
that contains `parent`, registered in `.GlobalEnv`. If `class` already
exists, stops with an error unless it already extends `parent`.

## Usage

``` r
defineOrCheckClass(class, parent)
```

## Arguments

- class:

  Character string — the S4 class name to define or check.

- parent:

  Character string — the required parent class.

## Value

Invisibly, `class`.
