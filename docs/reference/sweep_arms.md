# The order in which to project the scan values

With continuation each scan value starts from the attractor reached at
the previous one, so it pays to visit them in an order where consecutive
values are close together. When the value the model currently sits at is
known, that means working outwards from it in both directions; otherwise
the order the user gave is used, which is also what lets a decreasing
`scan_values` trace a hysteresis branch deliberately.

## Usage

``` r
sweep_arms(scan_values, current_scan_value = NULL)
```

## Arguments

- scan_values:

  The values to scan over.

- current_scan_value:

  The value the model currently sits at, or NULL.

## Value

A list of integer vectors, each holding indices into `scan_values` in
the order they should be projected. Each arm is to be started from the
unmodified model.

## Details

The two directions are returned as separate arms rather than as one
sequence, because each has to begin again at the model as it was given.
Run as one sequence they would carry the attractor from the far end of
the descending arm into the start of the ascending arm, which is the
opposite of starting each projection from a neighbour, and in a model
with coexisting attractors would follow the wrong branch.
