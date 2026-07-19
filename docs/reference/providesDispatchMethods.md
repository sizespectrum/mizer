# Does an installed extension register dispatch methods for its own class?

An extension package participates in dispatch by registering S3 methods
for mizer generics keyed on its marker class (e.g.
`getEncounter.mizerMR`). This checks the package namespace's own S3
method registry for any method whose class is the extension name or its
sim variant. Because S3 method registration does not require the S4
marker class to exist, this lets mizer recognise a dispatching extension
*before* creating its class, so extension packages no longer have to
define the marker class statically. Defining it statically as
`contains = "MizerParams"` would in fact prevent the package from being
chained with other extensions, since a sealed class cannot be
re-parented into the chain (see
[`defineExtensionClasses()`](https://sizespectrum.org/mizer/reference/defineExtensionClasses.md)).

## Usage

``` r
providesDispatchMethods(name)
```

## Arguments

- name:

  The extension identifier (its marker class / package name).

## Value

`TRUE` if the loaded namespace `name` registers S3 methods for class
`name` or `paste0(name, "Sim")`, otherwise `FALSE`.
