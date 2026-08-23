# Report the scan values that did not settle on a fixed point

Report the scan values that did not settle on a fixed point

## Usage

``` r
report_scan_convergence(
  scan_values,
  attractors,
  terminations,
  scan_name,
  t_max,
  t_sample
)
```

## Arguments

- scan_values:

  The values that were scanned.

- attractors:

  The attractor reached at each value, one per value.

- terminations:

  Why the run at each value stopped, one per value.

- scan_name:

  The name of the scanned quantity.

- t_max:

  The time limit that was used.

- t_sample:

  The averaging window that was used.

## Value

Nothing; called for its messages.
