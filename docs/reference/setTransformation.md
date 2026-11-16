# Register a joint model transformation

Local transformation declarations belong to
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
or
[`setProcess()`](https://sizespectrum.org/mizer/reference/setProcess.md).
Use a joint declaration when one hook must update data shared by several
owners. It replaces their local declarations for its operation;
overlapping owners or update paths are errors. Hooks run on the
constructed candidate before dependent rates are rebuilt or the complete
model is validated. Missing support is an error; `action = "none"`
explicitly declares that the core operation needs no additional update
from this provider.

## Usage

``` r
setTransformation(params, name, owners = NULL, operation = NULL, definition)
```

## Arguments

- params:

  A MizerParams object.

- name:

  Name identifying the joint declaration.

- owners:

  Canonical owners: consumers, resource, reproduction or component:name.

- operation:

  One of add_species, remove_species, reorder_species, rename_species,
  size_grid, scale_abundance or scale_rates.

- definition:

  A list with action (none or update), transform_fun (a named callback
  for update), paths (a list of literal character-vector paths) and
  serialisable settings. NULL removes the named joint declaration.

## Value

The updated MizerParams object.
