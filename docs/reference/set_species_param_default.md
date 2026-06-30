# Set a species parameter to a default value

If the species parameter does not yet exist in the species parameter
data frame, then create it and fill it with the default. Otherwise use
the default only to fill in any NAs. Optionally gives a message if the
parameter did not already exist. The signal has class
`info_about_default`.

## Usage

``` r
set_species_param_default(object, parname, default, message = NULL)
```

## Arguments

- object:

  Either a MizerParams object or a species parameter data frame

- parname:

  A string with the name of the species parameter to set

- default:

  A single default value or a vector with one default value for each
  species

- message:

  A string with a message to be issued when the parameter did not
  already exist

## Value

The `object` with an updated column in the species params data frame.
