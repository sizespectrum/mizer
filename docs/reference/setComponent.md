# Add a dynamical ecosystem component

By default, mizer models any number of size-resolved consumer species
and a single size-resolved resource spectrum. Your model may require
additional components, like for example detritus or carrion or multiple
resources or .... This function allows you to set up such components.

## Usage

``` r
setComponent(
  params,
  component,
  initial_value,
  dynamics_fun,
  encounter_fun,
  mort_fun,
  component_params,
  colour = "grey",
  linetype = "solid"
)

removeComponent(params, component)

getComponent(params, component)
```

## Arguments

- params:

  A MizerParams object

- component:

  Name of the component of interest. If missing, a list of all
  components will be returned.

- initial_value:

  Initial value of the component

- dynamics_fun:

  Name of function to calculate value at the next time step

- encounter_fun:

  Name of function to calculate contribution to encounter rate.
  Optional.

- mort_fun:

  Name of function to calculate contribution to the mortality rate.
  Optional.

- component_params:

  Object holding the parameters needed by the component functions. This
  could for example be a named list of parameters. Optional.

- colour:

  Line colour to use for the component in plots. Defaults to `"grey"`.

- linetype:

  Line type to use for the component in plots. Defaults to `"solid"`.

## Value

The updated MizerParams object

For `getComponent`: A list with the entries `initial_value`,
`dynamics_fun`, `encounter_fun`, `mort_fun`, `component_params` for the
requested component. If the requested component does not exist, `NULL`
is returned. If no `component` argument is given, then a list of lists
for all components is returned.

## Details

The component can be a number, a vector, an array, a list, or any other
data structure you like.

If you set a component with a new name, the new component will be added
to the existing components. If you set a component with an existing
name, the `initial_value` and `dynamics_fun` are overwritten, while the
optional `encounter_fun`, `mort_fun` and `component_params` are only
changed if the corresponding arguments are supplied. You can remove a
component with `removeComponent()`.

## See also

"Extending mizer":
[`vignette("extending-mizer", package = "mizer")`](https://sizespectrum.org/mizer/articles/extending-mizer.md)

Other extension tools:
[`NOther()`](https://sizespectrum.org/mizer/reference/NOther.md),
[`clearExtensionChain()`](https://sizespectrum.org/mizer/reference/clearExtensionChain.md),
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md),
[`getRegisteredExtensions()`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.md),
[`initialNOther<-()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md),
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md),
[`registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.md),
[`registerExtensions()`](https://sizespectrum.org/mizer/reference/registerExtensions.md),
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
