# Length of the zero-padded grid for the predation kernel convolutions

A fast Fourier transform evaluates a circular convolution, in which the
part of a kernel that falls off one end of the grid wraps around to the
other end.
[`convolve_pred_kernel()`](https://sizespectrum.org/mizer/reference/convolve_pred_kernel.md)
makes the convolution linear by padding the grid with zeros. The
predators are on the last `no_w` points of the `no_w_full` points of the
grid, so a predator and a prey larger than it are at most `no_w - 1`
points apart. In a circular convolution of length `no_k` such a pair is
wrapped onto offset `no_k - (k - n) >= no_k - no_w + 1`. If `m_max` is
the largest offset (from `0` to `no_w_full - 1`) at which the kernel
weights are non-zero, choosing `no_k >= max(no_w_full, no_w + m_max)`
places every wrapped pair beyond `m_max`, where the kernel weights are
zero, both for the encounter and diffusion rates (prey on the full grid,
predators on the last `no_w` points) and for the predation rate
(predators on the last `no_w` points, prey on the full grid). When
`m_max <= no_w_full - no_w` (the kernel of even the smallest predator
does not reach below `min_w_pp`), no padding beyond `no_w_full` is
needed. We round up to a length whose only prime factors are 2, 3 and 5,
for which the FFT is fast.

## Usage

``` r
pred_kernel_fft_length(no_w_full, no_w, edition, phis = NULL)
```

## Arguments

- no_w_full:

  The number of sizes in `w_full`.

- no_w:

  The number of sizes in `w`.

- edition:

  The defaults edition of the model, see
  [`defaults_edition()`](https://sizespectrum.org/mizer/reference/defaults_edition.md).

- phis:

  Optional matrix or list of matrices of kernel weights, each with
  `no_w_full` columns corresponding to offsets `0` to `no_w_full - 1`.
  If omitted, the kernel is assumed to be non-zero up to offset
  `no_w_full - 1`.

## Value

The length of the padded grid.

## Details

The padding came with defaults edition 3. A model of an earlier edition
is not padded: its length is `no_w_full`, so that its convolutions stay
circular and it keeps the results it had with earlier versions of mizer.
