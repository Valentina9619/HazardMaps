# Plot best-fitting marginal distribution by grid cell

Maps each cell by the marginal family selected by
\[hm_fit_marginals()\]. Reproduces Fig. 6 of Clavijo Mesa et al. (2026).

## Usage

``` r
hm_plot_best_dist(
  x,
  fits,
  bbox = NULL,
  title = NULL,
  show_map = TRUE,
  palette = NULL,
  coord = "fixed"
)
```

## Arguments

- x:

  An `hm_hazard` object after \[hm_standardize_coords()\].

- fits:

  The data.frame returned by \[hm_fit_marginals()\].

- bbox:

  Optional named numeric vector with \`lon_min\`, \`lon_max\`,
  \`lat_min\`, and \`lat_max\` to control the plot extent. If \`NULL\`
  (default), the extent is derived from the coordinate table with a
  small padding so edge points are not clipped.

- title:

  Optional plot title.

- show_map:

  Logical. If \`TRUE\`, adds a world map background.

- palette:

  Optional named character vector mapping distribution labels to
  colours. Names must match the values in \`fits\$best_dist_label\`
  (e.g. \`"Weibull"\`, \`"Gamma"\`, \`"Lognormal"\`, \`"Gumbel"\`,
  \`"Zero-trunc. Gaussian"\`, \`"Zero-trunc. Laplace"\`, \`"No fit"\`).
  If \`NULL\` (default), a built-in muted palette is used.

- coord:

  Coordinate display method. Either \`"fixed"\` or \`"quickmap"\`.

## Value

A \`ggplot\` object.

## See also

\[hm_fit_marginals()\], \[hm_plot_marginal_fit()\],
\[hm_plot_dist_frequency()\]

## Examples

``` r
if (FALSE) { # \dontrun{
x     <- hm_read_netcdf(f, var = "i10fg")
x     <- hm_standardize_coords(x)
x     <- hm_decode_time(x)
x_ext <- hm_select_extreme_events(x, threshold = 25)
X     <- hm_to_matrix(x_ext)
fits  <- hm_fit_marginals(X)

hm_plot_best_dist(x_ext, fits)
} # }
```
